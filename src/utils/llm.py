import os
import ast
import litellm
import time
from typing import Optional, Dict, Any, List
from dotenv import load_dotenv
from litellm import completion
from src.variables import OPENAI_MODELS, DEEPSEEK_MODELS
from src.variables import USER_PROMPT, SYS_PROMPT
from src.variables import USER_PROMPT_V4, SYS_PROMPT_V4
from src.variables import USER_PROMPT_OPENAI, SYS_PROMPT_OPENAI
from src.variables import USER_PROMPT_DEEPSEEK, SYS_PROMPT_DEEPSEEK
from src.variables import ADDON_PROMPT_7_MEMBER, USER_PROMPT_DEEPSEEK_V4
from src.variables import ERROR_MAP, PROTECTING_GROUP_CONTEXT
from src.cache import cache_results
from src.utils.utils_molecule import validity_check, detect_seven_member_rings
from src.utils.job_context import logger as context_logger
from src.utils.stability_checks import stability_checker
from src.utils.hallucination_checks import hallucination_checker
from src.protecting_group import mask_protecting_groups_multisymbol
from src.deprotecting_group import unmask_protecting_groups_multisymbol, get_protecting_group_info
from src.utils.tools import get_available_tools, execute_tool

# Langfuse SDK for direct integration
try:
    from langfuse import Langfuse
    langfuse_client = Langfuse()
    LANGFUSE_AVAILABLE = True
except ImportError:
    LANGFUSE_AVAILABLE = False
    langfuse_client = None

load_dotenv()

# set the success callback to langfuse for logging
litellm.success_callback = ["langfuse"]
litellm.drop_params = True

metadata = {
    "generation_name": "prod",  # set langfuse generation name
    "project": "Retrosynthesis",  # set langfuse project name
    "version": "0.0.3",  # set langfuse version
    "trace_name": "prod",  # set langfuse Trace Name
    "trace_user_id": "sv",  # set langfuse Trace User ID
    "session_id": "prod",  # set langfuse Session ID
}
ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING",
                                    "true").lower() == "false" else True


def log_message(message: str, logger=None):
    """Log the message

    Parameters
    ----------
    message : str
        The message to be logged
    logger : _type_, optional
        The logger object, by default None

    Returns
    -------
    None
    """
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def obtain_prompt(LLM: str):
    """Obtain the prompt based on the LLM model

    Parameters
    ----------
    LLM : str
        The LLM model to be used

    Returns
    -------
    str, str, int
        The system prompt, user prompt and max completion tokens
    """
    advanced_prompt = False
    detector = LLM.split(":")
    if len(detector) > 1 and detector[1] == "adv":
        advanced_prompt = True
    print(f"Advanced Prompt: {advanced_prompt}")
    if advanced_prompt:
        if LLM in DEEPSEEK_MODELS:
            sys_prompt_final = SYS_PROMPT_V4
            user_prompt_final = USER_PROMPT_DEEPSEEK_V4
            max_completion_tokens = 8192 * 2
        elif LLM in OPENAI_MODELS:
            sys_prompt_final = SYS_PROMPT_OPENAI
            user_prompt_final = USER_PROMPT_OPENAI
            max_completion_tokens = 8192
        else:
            sys_prompt_final = SYS_PROMPT_V4
            user_prompt_final = USER_PROMPT_V4
            max_completion_tokens = 4096
    else:
        if LLM in DEEPSEEK_MODELS:
            sys_prompt_final = SYS_PROMPT_DEEPSEEK
            user_prompt_final = USER_PROMPT_DEEPSEEK
            max_completion_tokens = 8192 * 2
        elif LLM in OPENAI_MODELS:
            sys_prompt_final = SYS_PROMPT_OPENAI
            user_prompt_final = USER_PROMPT_OPENAI
            max_completion_tokens = 8192
        else:
            sys_prompt_final = SYS_PROMPT
            user_prompt_final = USER_PROMPT
            max_completion_tokens = 4096
    return sys_prompt_final, user_prompt_final, max_completion_tokens


@cache_results
def call_LLM(molecule: str,
             LLM: str = "claude-opus-4-20250514",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None,
             use_protecting_group_feature: bool = False) -> tuple[int, str]:
    """Calls the LLM model to predict the next step

    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        The LLM model to be used, by default "claude-opus-4-20250514"
    temperature : float, optional
        The temperature for sampling, by default 0.0
    messages : Optional[list[dict]], optional
        The conversation history, by default None

    Returns
    -------
    tuple[int, str]
        The status code and the response text
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    log_message(f"Calling {LLM} with molecule: {molecule}", logger)

    # Check for seven-member rings
    if detect_seven_member_rings(molecule):
        log_message(f"Detected seven member ring in molecule: {molecule}",
                    logger)
        add_on = ADDON_PROMPT_7_MEMBER
    else:
        add_on = ""

    # Check for protecting groups and add context
    if use_protecting_group_feature:
        masked_smiles = mask_protecting_groups_multisymbol(molecule)
        if masked_smiles != molecule and masked_smiles != "INVALID_SMILES":
            log_message(
                f"Detected protecting groups in molecule: {molecule} -> {masked_smiles}",
                logger)
            protecting_group_context = PROTECTING_GROUP_CONTEXT.format(
                molecule=molecule, masked_smiles=masked_smiles)
            add_on += protecting_group_context

    sys_prompt_final, user_prompt_final, max_completion_tokens = obtain_prompt(
        LLM)
    LLM = LLM.split(":")[0]

    params = {
        "model": LLM,
        "max_completion_tokens": max_completion_tokens,
        "temperature": temperature,
        "seed": 42,
        "top_p": 0.9,
        "metadata": metadata,
    }

    if LLM in DEEPSEEK_MODELS:
        user_prompt_final += add_on

    if "3-7" in LLM:
        params["max_tokens"] = 13192 + 5000
        params["temperature"] = 1
        params.pop("top_p", None)
        params.pop("max_completion_tokens", None)
        params['thinking'] = {"type": "enabled", "budget_tokens": 5000}

    if messages is None:
        messages = [{
            "role": "system",
            "content": sys_prompt_final + add_on
        }, {
            "role":
            "user",
            "content":
            user_prompt_final.replace('{target_smiles}', molecule) + "\n\n" +
            add_on
        }]
    params["messages"] = messages

    try:
        # Call the LLM model
        response = completion(**params)

        res_text = response.choices[0].message.content
    except Exception as e:
        log_message(f"Error in calling {LLM}: {e}", logger)
        log_message(f"Retrying call to {LLM}", logger)
        try:
            response = completion(**params)
            res_text = response.choices[0].message.content
        except Exception as e:
            log_message(f"2nd Error in calling {LLM}: {e}", logger)
            log_message(f"Exiting call to {LLM}", logger)
            return 400, ""
    log_message(f"Received response from LLM: {res_text}", logger)
    return 200, res_text


@cache_results
def agentic_call_LLM(
        molecule: str,
        LLM: str = "claude-opus-4-20250514",
        temperature: float = 0.0,
        messages: Optional[list[dict]] = None,
        use_protecting_group_feature: bool = False,
        max_iterations: int = 5,
        trace_id: Optional[str] = None) -> tuple[int, str, Dict[str, Any]]:
    """Agentic LLM call with tool use support and Langfuse logging
    
    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        The LLM model to be used, by default "claude-opus-4-20250514"
    temperature : float, optional
        The temperature for sampling, by default 0.0
    messages : Optional[list[dict]], optional
        The conversation history, by default None
    use_protecting_group_feature : bool, optional
        Whether to use protecting group feature, by default False
    max_iterations : int, optional
        Maximum number of agent iterations, by default 5
    trace_id : Optional[str], optional
        Langfuse trace ID for hierarchical logging, by default None
    
    Returns
    -------
    tuple[int, str, Dict[str, Any]]
        The status code, response text, and tool usage statistics
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    log_message(f"Calling agentic {LLM} with molecule: {molecule}", logger)

    # Initialize tool usage statistics
    tool_stats = {
        "total_iterations": 0,
        "total_tool_calls": 0,
        "tool_calls_by_name": {},
        "successful_validations": 0,
        "failed_validations": 0
    }

    # Initialize Langfuse trace
    trace = None
    if LANGFUSE_AVAILABLE and langfuse_client:
        try:
            trace = langfuse_client.trace(id=trace_id,
                                          name="agentic_retrosynthesis",
                                          user_id="sv",
                                          metadata={
                                              "molecule":
                                              molecule,
                                              "model":
                                              LLM,
                                              "temperature":
                                              temperature,
                                              "max_iterations":
                                              max_iterations,
                                              "use_protecting_groups":
                                              use_protecting_group_feature
                                          })
        except Exception as e:
            log_message(f"Langfuse trace creation failed: {e}", logger)

    # Check for seven-member rings
    if detect_seven_member_rings(molecule):
        log_message(f"Detected seven member ring in molecule: {molecule}",
                    logger)
        add_on = ADDON_PROMPT_7_MEMBER
    else:
        add_on = ""

    # Check for protecting groups and add context
    if use_protecting_group_feature:
        masked_smiles = mask_protecting_groups_multisymbol(molecule)
        if masked_smiles != molecule and masked_smiles != "INVALID_SMILES":
            log_message(
                f"Detected protecting groups in molecule: {molecule} -> {masked_smiles}",
                logger)
            protecting_group_context = PROTECTING_GROUP_CONTEXT.format(
                molecule=molecule, masked_smiles=masked_smiles)
            add_on += protecting_group_context

    sys_prompt_final, user_prompt_final, max_completion_tokens = obtain_prompt(
        LLM)
    LLM = LLM.split(":")[0]

    params = {
        "model": LLM,
        "max_completion_tokens": max_completion_tokens,
        "temperature": temperature,
        "seed": 42,
        "top_p": 0.9,
        "metadata": metadata,
    }

    if LLM in DEEPSEEK_MODELS:
        user_prompt_final += add_on

    if "3-7" in LLM:
        params["max_tokens"] = 13192 + 5000
        params["temperature"] = 1
        params.pop("top_p", None)
        params.pop("max_completion_tokens", None)
        params['thinking'] = {"type": "enabled", "budget_tokens": 5000}

    # Prepare initial messages
    if messages is None:
        messages = [{
            "role": "system",
            "content": sys_prompt_final + add_on
        }, {
            "role":
            "user",
            "content":
            user_prompt_final.replace('{target_smiles}', molecule) + "\n\n" +
            add_on
        }]

    # Get available tools
    tools = get_available_tools()

    # Agent loop
    iteration = 0
    conversation_messages = messages.copy()

    while iteration < max_iterations:
        tool_stats["total_iterations"] = iteration + 1
        iteration_start_time = time.time()

        # Create iteration span in Langfuse
        iteration_span = None
        if trace:
            try:
                iteration_span = trace.span(
                    name=f"agent_iteration_{iteration + 1}",
                    metadata={
                        "iteration": iteration + 1,
                        "max_iterations": max_iterations
                    })
            except Exception as e:
                log_message(f"Langfuse iteration span creation failed: {e}",
                            logger)

        log_message(f"Agent iteration {iteration + 1}/{max_iterations}",
                    logger)

        # Prepare params for this iteration
        iter_params = params.copy()
        iter_params["messages"] = conversation_messages
        iter_params["tools"] = tools

        try:
            # Create Langfuse generation
            generation = None
            if iteration_span:
                try:
                    generation = iteration_span.generation(
                        name="llm_call",
                        model=LLM,
                        input=conversation_messages,
                        metadata=params)
                except Exception as e:
                    log_message(f"Langfuse generation creation failed: {e}",
                                logger)

            # Call the LLM
            response = completion(**iter_params)

            # Log generation output
            if generation:
                try:
                    generation.end(
                        output=response.choices[0].message.content if hasattr(
                            response.choices[0].message, 'content') else None,
                        usage={
                            "input_tokens":
                            getattr(response.usage, 'input_tokens', 0)
                            if hasattr(response, 'usage') else 0,
                            "output_tokens":
                            getattr(response.usage, 'output_tokens', 0)
                            if hasattr(response, 'usage') else 0,
                            "total_tokens":
                            getattr(response.usage, 'total_tokens', 0)
                            if hasattr(response, 'usage') else 0
                        },
                        metadata={
                            "stop_reason":
                            response.choices[0].finish_reason
                            if response.choices else "unknown"
                        })
                except Exception as e:
                    log_message(f"Langfuse generation end failed: {e}", logger)

            res_message = response.choices[0].message
            stop_reason = response.choices[0].finish_reason

            log_message(f"LLM response stop_reason: {stop_reason}", logger)

            # Add assistant response to conversation
            # Include tool_calls if present for OpenAI format
            assistant_msg = {
                "role":
                "assistant",
                "content":
                res_message.content if hasattr(res_message, 'content')
                and res_message.content else None
            }
            if hasattr(res_message, 'tool_calls') and res_message.tool_calls:
                assistant_msg["tool_calls"] = res_message.tool_calls

            conversation_messages.append(assistant_msg)

            # Check if there are tool calls
            if stop_reason == "tool_use" or (hasattr(res_message, 'tool_calls')
                                             and res_message.tool_calls):
                # Handle tool calls
                tool_calls = []

                # Extract tool calls (format varies by provider)
                if hasattr(res_message,
                           'tool_calls') and res_message.tool_calls:
                    tool_calls = res_message.tool_calls

                log_message(f"Processing {len(tool_calls)} tool calls", logger)

                tool_results = []
                for tool_call in tool_calls:
                    tool_stats["total_tool_calls"] += 1

                    # Extract tool info
                    tool_name = tool_call.function.name if hasattr(
                        tool_call, 'function') else getattr(
                            tool_call, 'name', 'unknown')
                    tool_input = tool_call.function.arguments if hasattr(
                        tool_call, 'function') else getattr(
                            tool_call, 'input', {})
                    tool_id = getattr(tool_call, 'id',
                                      f"call_{iteration}_{len(tool_results)}")

                    # Parse arguments if string
                    if isinstance(tool_input, str):
                        try:
                            tool_input = ast.literal_eval(tool_input)
                        except:
                            try:
                                import json
                                tool_input = json.loads(tool_input)
                            except:
                                tool_input = {"raw": tool_input}

                    # Track tool usage
                    tool_stats["tool_calls_by_name"][
                        tool_name] = tool_stats["tool_calls_by_name"].get(
                            tool_name, 0) + 1

                    log_message(
                        f"Executing tool: {tool_name} with input: {tool_input}",
                        logger)

                    # Create tool span in Langfuse
                    tool_span = None
                    if iteration_span:
                        try:
                            tool_span = iteration_span.span(
                                name=f"tool_call_{tool_name}",
                                input=tool_input,
                                metadata={
                                    "tool_name": tool_name,
                                    "tool_id": tool_id
                                })
                        except Exception as e:
                            log_message(
                                f"Langfuse tool span creation failed: {e}",
                                logger)

                    # Execute tool
                    tool_result = execute_tool(tool_name, tool_input)

                    # Track validation results
                    if tool_name == "validate_smiles":
                        if tool_result.get("valid", False):
                            tool_stats["successful_validations"] += 1
                        else:
                            tool_stats["failed_validations"] += 1

                    log_message(f"Tool result: {tool_result}", logger)

                    # End tool span
                    if tool_span:
                        try:
                            tool_span.end(output=tool_result,
                                          metadata={
                                              "success":
                                              tool_result.get("success", True),
                                              "error":
                                              tool_result.get("error", None)
                                          })
                        except Exception as e:
                            log_message(f"Langfuse tool span end failed: {e}",
                                        logger)

                    # Format for OpenAI/litellm compatibility
                    tool_results.append({
                        "role": "tool",
                        "tool_call_id": tool_id,
                        "content": str(tool_result)
                    })

                # Add tool results to conversation
                # For OpenAI format, add each tool result as separate message
                for tool_result_msg in tool_results:
                    conversation_messages.append(tool_result_msg)

                # End iteration span
                if iteration_span:
                    try:
                        iteration_span.end(
                            metadata={
                                "tool_calls":
                                len(tool_calls),
                                "iteration_time":
                                time.time() - iteration_start_time
                            })
                    except Exception as e:
                        log_message(f"Langfuse iteration span end failed: {e}",
                                    logger)

                # Continue to next iteration
                iteration += 1
                continue

            # No tool use, we have final answer
            res_text = res_message.content if hasattr(res_message,
                                                      'content') else ""

            # End iteration span
            if iteration_span:
                try:
                    iteration_span.end(
                        metadata={
                            "final_iteration": True,
                            "iteration_time": time.time() -
                            iteration_start_time
                        })
                except Exception as e:
                    log_message(f"Langfuse iteration span end failed: {e}",
                                logger)

            # Log final result event
            if trace:
                try:
                    trace.event(name="agent_completed",
                                metadata={
                                    "total_iterations":
                                    tool_stats["total_iterations"],
                                    "total_tool_calls":
                                    tool_stats["total_tool_calls"],
                                    "successful_validations":
                                    tool_stats["successful_validations"],
                                    "failed_validations":
                                    tool_stats["failed_validations"]
                                })
                except Exception as e:
                    log_message(f"Langfuse event creation failed: {e}", logger)

            log_message(f"Agent completed in {iteration + 1} iterations",
                        logger)
            log_message(f"Tool usage statistics: {tool_stats}", logger)
            log_message(f"Received final response from LLM: {res_text}",
                        logger)

            return 200, res_text, tool_stats

        except Exception as e:
            log_message(
                f"Error in agentic LLM call iteration {iteration + 1}: {e}",
                logger)

            # End iteration span with error
            if iteration_span:
                try:
                    iteration_span.end(
                        metadata={
                            "error": str(e),
                            "iteration_time": time.time() -
                            iteration_start_time
                        })
                except Exception as ex:
                    log_message(f"Langfuse iteration span end failed: {ex}",
                                logger)

            if iteration == 0:
                # First iteration failed, retry once
                log_message(f"Retrying agentic call to {LLM}", logger)
                try:
                    # Simplified retry without tools
                    retry_params = params.copy()
                    retry_params["messages"] = messages
                    response = completion(**retry_params)
                    res_text = response.choices[0].message.content
                    return 200, res_text, tool_stats
                except Exception as e2:
                    log_message(f"2nd Error in agentic call to {LLM}: {e2}",
                                logger)
                    return 400, "", tool_stats
            else:
                return 400, "", tool_stats

    # Max iterations reached
    log_message(
        f"Max iterations ({max_iterations}) reached without final answer",
        logger)

    if trace:
        try:
            trace.event(name="max_iterations_reached",
                        metadata={
                            "total_iterations": tool_stats["total_iterations"],
                            "total_tool_calls": tool_stats["total_tool_calls"]
                        })
        except Exception as e:
            log_message(f"Langfuse event creation failed: {e}", logger)

    # Try to extract last assistant message
    for msg in reversed(conversation_messages):
        if msg.get("role") == "assistant" and msg.get("content"):
            return 200, msg["content"], tool_stats

    return 400, "Max iterations reached without final answer", tool_stats


def split_cot_json(res_text: str) -> tuple[int, list[str], str]:
    """Parse the LLM response to extract the thinking steps and json content

    Parameters
    ----------
    res_text : str
        The response text from the LLM model

    Returns
    -------
    tuple[int, list[str], str]
        The status code, thinking steps and json content
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        # extract the content within <cot> </cot> tags as thinking content
        thinking_content = res_text[res_text.find("<cot>\n") +
                                    6:res_text.find("</cot>")]
        if not thinking_content:
            return 501, [], ""

        # split the thinking content into individual steps based on the <thinking> </thinking> tags
        thinking_steps = thinking_content.split("<thinking")[1:]
        thinking_steps = [
            step[:step.find("</thinking>")] for step in thinking_steps
        ]
        if not thinking_steps:
            return 501, [], ""
    except Exception as e:
        log_message(f"Error in parsing obtaining COT: {e}", logger)
        return 501, [], ""

    try:
        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
        if not json_content:
            return 501, [], ""
    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 501, [], ""
    return 200, thinking_steps, json_content


def split_json_openAI(res_text: str) -> tuple[int, str]:
    """Split the response text from OpenAI models to extract the molecules
    Note: OpenAI O-series models do not provide Chain of Thoughts (COT) in the response

    Parameters
    ----------
    res_text : str
        The response text from the OpenAI model

    Returns
    -------
    tuple[int, str]
        the status code and json content
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
        if not json_content:
            return 502, ""

    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 502, ""
    return 200, json_content


def split_json_deepseek(res_text: str) -> tuple[int, list[str], str]:
    """Parse the LLM response to extract the thinking steps and json content

    Parameters
    ----------
    res_text : str
        The response text from the LLM model

    Returns
    -------
    tuple[int, list[str], str]
        The status code, thinking steps and json content
    """
    logger = context_logger.get() if ENABLE_LOGGING else None

    try:
        # extract the content within <cot> </cot> tags as thinking content
        thinking_content = res_text[res_text.find("<think>\n") +
                                    6:res_text.find("</think>")]
        if not thinking_content:
            return 503, [], ""

        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
        if not json_content:
            return 503, [], ""

    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 503, [], ""
    return 200, [thinking_content], json_content


def split_json_master(res_text: str, model: str) -> tuple[int, list[str], str]:
    """Split the response text based on the model

    Parameters
    ----------
    res_text : str
        The response text from the LLM model
    model : str
        The LLM model used

    Returns
    -------
    tuple[int, list[str], str]
        The status code, thinking steps and json content
    """
    try:
        if model in DEEPSEEK_MODELS:
            status_code, thinking_steps, json_content = split_json_deepseek(
                res_text)
        elif model in OPENAI_MODELS:
            status_code, json_content = split_json_openAI(res_text)
            thinking_steps = []
        else:
            status_code, thinking_steps, json_content = split_cot_json(
                res_text)
    except Exception as e:
        return 505, [], ""

    return status_code, thinking_steps, json_content


def validate_split_json(
        json_content: str) -> tuple[int, list[str], list[str], list[int]]:
    """Validate the split json content from LLM response

    Parameters
    ----------
    json_content : str
        The json content from the LLM response

    Returns
    -------
    tuple[int, list[str], list[str], list[int]]
        The status code, list of molecules, list of explanations and list of confidence scores
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        result_list = ast.literal_eval(json_content)
        res_molecules = result_list['data']
        res_explanations = result_list['explanation']
        res_confidence = result_list['confidence_scores']
    except Exception as e:
        log_message(f"Error in parsing response: {e}", logger)
        return 504, [], [], []
    return 200, res_molecules, res_explanations, res_confidence


def llm_pipeline(
    molecule: str,
    LLM: str = "claude-opus-4-20250514",
    messages: Optional[list[dict]] = None,
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False
) -> tuple[list[list[str]], list[str], list[float]]:
    """Pipeline to call LLM and validate the results

    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        LLM to be used for retrosynthesis , by default "claude-opus-4-20250514"
    messages : Optional[list[dict]], optional
        Conversation history, by default None

    Returns
    -------
    tuple[list[list[str]], list[str], list[float]]
        The output pathways, explanations and confidence scores
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    output_pathways: list[list[str]] = []
    output_explanations: list[str] = []
    output_confidence: list[float] = []
    run = 0.0
    if stability_flag.lower() == "true" or hallucination_check.lower(
    ) == "true":
        max_run = 1.5
    else:
        max_run = 0.6
    while (output_pathways == [] and run < max_run):
        log_message(f"Calling LLM with molecule: {molecule} and run: {run}",
                    logger)

        # Selecting the model based on the run number
        current_model = LLM
        if LLM in DEEPSEEK_MODELS and run > 0.0:
            current_model = "claude-opus-4-20250514"

        # --------------------
        # Call LLM
        status_code, res_text = call_LLM(
            molecule,
            current_model,
            messages=messages,
            temperature=run,
            use_protecting_group_feature=use_protecting_group_feature)
        if status_code != 200:
            log_message(f"Error in calling LLM: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Split the response text
        status_code, thinking_steps, json_content = split_json_master(
            res_text, current_model)
        if status_code != 200:
            log_message(f"Error in splitting cot json: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Validate the split json content
        status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
            json_content)
        if status_code != 200:
            log_message(f"Error in validating split json content: {res_text}",
                        logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Check the validity of the molecules obtained from LLM
        output_pathways, output_explanations, output_confidence = validity_check(
            molecule, res_molecules, res_explanations, res_confidence)

        # --------------------
        # Stability check
        if stability_flag.lower() == "true":
            status_code, stable_pathways = stability_checker(output_pathways)
            if status_code != 200:
                log_message(f"Error in stability check: {stable_pathways}",
                            logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = stable_pathways

        # --------------------
        # Hallucination check
        if hallucination_check.lower() == "true":
            log_message(
                f"Calling hallucination check with pathways: {output_pathways}",
                logger)
            status_code, hallucination_pathways = hallucination_checker(
                molecule, output_pathways)
            if status_code != 200:
                log_message(
                    f"Error in hallucination check: {hallucination_pathways}",
                    logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = hallucination_pathways

        # --------------------
        # Update the messages for the next call
        # if output_pathways == []:
        #     messages = [{
        #         "role": "system",
        #         "content": SYS_PROMPT
        #     }, {
        #         "role": "user",
        #         "content": USER_PROMPT
        #     }, {
        #         "role": "assistant",
        #         "content": res_text
        #     }, {
        #         "role": "user",
        #         "content": "<add something here>"
        #     }]

        log_message(
            f"Output Pathways: {output_pathways},\n\
                Output Explanations: {output_explanations},\n\
                    Output Confidence: {output_confidence}", logger)
        run += 0.1

    return output_pathways, output_explanations, output_confidence


def agentic_llm_pipeline(
    molecule: str,
    LLM: str = "claude-opus-4-20250514",
    messages: Optional[list[dict]] = None,
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
    max_tool_iterations: int = 5,
    trace_id: Optional[str] = None
) -> tuple[list[list[str]], list[str], list[float], Dict[str, Any]]:
    """Agentic pipeline to call LLM with tool support and validate the results
    
    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        LLM to be used for retrosynthesis, by default "claude-opus-4-20250514"
    messages : Optional[list[dict]], optional
        Conversation history, by default None
    stability_flag : str, optional
        Whether to run stability checks, by default "False"
    hallucination_check : str, optional
        Whether to run hallucination checks, by default "False"
    use_protecting_group_feature : bool, optional
        Whether to use protecting group feature, by default False
    max_tool_iterations : int, optional
        Maximum number of agent iterations, by default 5
    trace_id : Optional[str], optional
        Langfuse trace ID for hierarchical logging, by default None
    
    Returns
    -------
    tuple[list[list[str]], list[str], list[float], Dict[str, Any]]
        The output pathways, explanations, confidence scores, and tool usage statistics
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    output_pathways: list[list[str]] = []
    output_explanations: list[str] = []
    output_confidence: list[float] = []
    all_tool_stats: Dict[str, Any] = {
        "runs": [],
        "total_tool_calls": 0,
        "total_iterations": 0,
        "successful_validations": 0,
        "failed_validations": 0
    }

    run = 0.0
    if stability_flag.lower() == "true" or hallucination_check.lower(
    ) == "true":
        max_run = 1.5
    else:
        max_run = 0.6

    while (output_pathways == [] and run < max_run):
        log_message(
            f"Calling agentic LLM with molecule: {molecule} and run: {run}",
            logger)

        # Selecting the model based on the run number
        current_model = LLM
        if LLM in DEEPSEEK_MODELS and run > 0.0:
            current_model = "claude-opus-4-20250514"

        # --------------------
        # Call Agentic LLM
        status_code, res_text, tool_stats = agentic_call_LLM(
            molecule,
            current_model,
            messages=messages,
            temperature=run,
            use_protecting_group_feature=use_protecting_group_feature,
            max_iterations=max_tool_iterations,
            trace_id=trace_id)

        # Track tool statistics
        all_tool_stats["runs"].append({
            "run": run,
            "temperature": run,
            "tool_stats": tool_stats
        })
        all_tool_stats["total_tool_calls"] += tool_stats.get(
            "total_tool_calls", 0)
        all_tool_stats["total_iterations"] += tool_stats.get(
            "total_iterations", 0)
        all_tool_stats["successful_validations"] += tool_stats.get(
            "successful_validations", 0)
        all_tool_stats["failed_validations"] += tool_stats.get(
            "failed_validations", 0)

        if status_code != 200:
            log_message(f"Error in calling agentic LLM: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Split the response text
        status_code, thinking_steps, json_content = split_json_master(
            res_text, current_model)
        if status_code != 200:
            log_message(f"Error in splitting cot json: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Validate the split json content
        status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
            json_content)
        if status_code != 200:
            log_message(f"Error in validating split json content: {res_text}",
                        logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Check the validity of the molecules obtained from LLM
        output_pathways, output_explanations, output_confidence = validity_check(
            molecule, res_molecules, res_explanations, res_confidence)

        # --------------------
        # Stability check
        if stability_flag.lower() == "true":
            status_code, stable_pathways = stability_checker(output_pathways)
            if status_code != 200:
                log_message(f"Error in stability check: {stable_pathways}",
                            logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = stable_pathways

        # --------------------
        # Hallucination check
        if hallucination_check.lower() == "true":
            log_message(
                f"Calling hallucination check with pathways: {output_pathways}",
                logger)
            status_code, hallucination_pathways = hallucination_checker(
                molecule, output_pathways)
            if status_code != 200:
                log_message(
                    f"Error in hallucination check: {hallucination_pathways}",
                    logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = hallucination_pathways

        log_message(
            f"Output Pathways: {output_pathways},\n\
                Output Explanations: {output_explanations},\n\
                    Output Confidence: {output_confidence}", logger)
        log_message(f"Tool usage statistics: {all_tool_stats}", logger)
        run += 0.1

    return output_pathways, output_explanations, output_confidence, all_tool_stats


def get_error_log(status_code: int):
    """Prints error message based on the status code.

    Parameters
    ----------
    status_code : int
        Status Code of the error.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None

    if status_code in ERROR_MAP:
        description = ERROR_MAP[status_code]
        log_message(f"Error Code: {status_code},\n Description: {description}",
                    logger)
    else:
        log_message(f"Error Code: {status_code} is not recognized.", logger)
