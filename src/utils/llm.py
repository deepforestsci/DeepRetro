import os
import ast
import litellm
from typing import Optional
from dotenv import load_dotenv
from litellm import completion
from src.variables import OPENAI_MODELS, DEEPSEEK_MODELS
from src.variables import USER_PROMPT, SYS_PROMPT
from src.variables import USER_PROMPT_V4, SYS_PROMPT_V4
from src.variables import USER_PROMPT_OPENAI, SYS_PROMPT_OPENAI
from src.variables import USER_PROMPT_DEEPSEEK, SYS_PROMPT_DEEPSEEK
from src.variables import ADDON_PROMPT_7_MEMBER, USER_PROMPT_DEEPSEEK_V4
from src.variables import ERROR_MAP
from src.cache import cache_results
from src.utils.utils_molecule import validity_check, detect_seven_member_rings
from src.utils.job_context import logger as context_logger
from src.utils.stability_checks import stability_checker
from src.utils.hallucination_checks import hallucination_checker

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
             LLM: str = "claude-3-opus-20240229",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None) -> tuple[int, str]:
    """Calls the LLM model to predict the next step

    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        The LLM model to be used, by default "claude-3-opus-20240229"
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

    if detect_seven_member_rings(molecule):
        log_message(f"Detected seven member ring in molecule: {molecule}",
                    logger)
        add_on = ADDON_PROMPT_7_MEMBER
    else:
        add_on = ""

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
    LLM: str = "claude-3-opus-20240229",
    messages: Optional[list[dict]] = None,
    stability_flag: str = "False",
    hallucination_check: str = "False",
    feedback_flag: str = "False"
) -> tuple[list[list[str]], list[str], list[float]]:
    """Pipeline to call LLM and validate the results

    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        LLM to be used for retrosynthesis , by default "claude-3-opus-20240229"
    messages : Optional[list[dict]], optional
        Conversation history, by default None
    stability_flag : str, optional
        Flag to enable stability check, by default "False"
    hallucination_check : str, optional
        Flag to enable hallucination check, by default "False"
    feedback_flag : str, optional
        Flag to enable feedback mechanism, by default "False"

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
        max_run = 1.2
    else:
        max_run = 0.6

    # Store the last response text for feedback
    last_res_text = ""

    while (output_pathways == [] and run < max_run):
        log_message(f"Calling LLM with molecule: {molecule} and run: {run}",
                    logger)

        # Selecting the model based on the run number
        current_model = LLM
        if LLM in DEEPSEEK_MODELS and run > 0.0:
            current_model = "anthropic/claude-3-7-sonnet-20250219"

        # If feedback is enabled and this is not the first run, provide feedback
        if feedback_flag.lower(
        ) == "true" and run > 0.0 and messages is not None and last_res_text:
            feedback_message = generate_feedback(last_res_text, current_model,
                                                 molecule)
            if feedback_message:
                # Update messages with feedback for the next call
                messages = update_messages_with_feedback(
                    messages, feedback_message, last_res_text)
                log_message(f"Added feedback: {feedback_message}", logger)

        # --------------------
        # Call LLM
        if run == 0.0 or feedback_flag.lower() == "false":
            messages = None

        status_code, res_text = call_LLM(molecule=molecule,
                                         LLM=current_model,
                                         messages=messages,
                                         temperature=run)

        # Store the response for potential feedback
        last_res_text = res_text

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

        validity_failed = (output_pathways == [])
        stability_failed = False
        hallucination_failed = False

        # --------------------
        # Stability check
        if stability_flag.lower() == "true" and output_pathways != []:
            status_code, stable_pathways = stability_checker(output_pathways)
            if status_code != 200:
                log_message(f"Error in stability check: {stable_pathways}",
                            logger)
                run += 0.1
                get_error_log(status_code)
                continue

            # Check if stability filtering removed all pathways
            stability_failed = (stable_pathways == []
                                and output_pathways != [])
            output_pathways = stable_pathways

        # --------------------
        # Hallucination check
        if hallucination_check.lower() == "true" and output_pathways != []:
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

            # Check if hallucination filtering removed all pathways
            hallucination_failed = (hallucination_pathways == []
                                    and output_pathways != [])
            output_pathways = hallucination_pathways

        # If we're going to do another run and feedback is enabled, prepare messages
        if output_pathways == [] and feedback_flag.lower() == "true":
            # Initialize messages if not provided
            if messages is None:
                sys_prompt_final, user_prompt_final, _ = obtain_prompt(
                    current_model)

                # Check for seven member rings
                add_on = ""
                if detect_seven_member_rings(molecule):
                    add_on = ADDON_PROMPT_7_MEMBER

                messages = [{
                    "role": "system",
                    "content": sys_prompt_final + add_on
                }, {
                    "role":
                    "user",
                    "content":
                    user_prompt_final.replace('{target_smiles}', molecule) +
                    "\n\n" + add_on
                }, {
                    "role": "assistant",
                    "content": res_text
                }]

        log_message(
            f"Output Pathways: {output_pathways},\n\
                Output Explanations: {output_explanations},\n\
                    Output Confidence: {output_confidence}", logger)
        run += 0.1

    return output_pathways, output_explanations, output_confidence


def generate_feedback(response_text: str, model: str, molecule: str) -> str:
    """
    Generate feedback based on the most recent failed response.
    
    Parameters
    ----------
    response_text : str
        The LLM response text that failed validation
    model : str
        The model that generated the response
    molecule : str
        The target molecule
        
    Returns
    -------
    str
        Specific feedback about what went wrong
    """
    feedback = "Your previous attempt had issues. "

    # Parse the response to identify issues
    try:
        # Extract JSON content
        if model in DEEPSEEK_MODELS:
            json_content = response_text[response_text.find("<json>\n") +
                                         7:response_text.find("</json>")]
        elif model in OPENAI_MODELS:
            json_content = response_text[response_text.find("<json>\n") +
                                         7:response_text.find("</json>")]
        else:  # Claude models
            json_content = response_text[response_text.find("<json>\n") +
                                         7:response_text.find("</json>")]

        # Parse the JSON content
        try:
            result_list = ast.literal_eval(json_content)
            res_molecules = result_list.get('data', [])

            # Check if molecules are valid
            for idx, pathway in enumerate(res_molecules):
                invalid_smiles = []
                for smiles in pathway:
                    # Use simple validity check - more detailed checking happens in validity_check function
                    # This is just for feedback purposes
                    if not smiles or '.' in smiles or '*' in smiles:
                        invalid_smiles.append(smiles)

                if invalid_smiles:
                    feedback += f"Pathway {idx+1} contained invalid SMILES: {invalid_smiles}. "

            # If no specific issues found but still failed, provide general guidance
            if "invalid" not in feedback:
                feedback += (
                    "Please check that your proposed reactants are commercially available or easily synthesizable. "
                    "Ensure all molecules have valid structures with proper valence. "
                )
        except (SyntaxError, ValueError):
            feedback += "Your JSON format was incorrect or could not be parsed. Please ensure you provide valid JSON in the <json> tags. "
    except Exception:
        feedback += "Please ensure your response includes properly formatted <json> and </json> tags with valid content. "

    feedback += (
        f"Recall that we're trying to perform retrosynthesis on this molecule: {molecule}. "
        "Focus on chemically feasible transformations, valid molecular structures, and "
        "commercially available starting materials. ")

    return feedback


def update_messages_with_feedback(messages: list[dict], feedback: str,
                                  last_response: str) -> list[dict]:
    """
    Update the message history with feedback for the next LLM call.
    
    Parameters
    ----------
    messages : list[dict]
        The current message history
    feedback : str
        The feedback to add
    last_response : str
        The last response from the LLM
        
    Returns
    -------
    list[dict]
        Updated message history with feedback
    """
    # If the last message was from the assistant, add the user feedback
    if messages and messages[-1]["role"] == "assistant":
        messages.append({"role": "user", "content": feedback})
    # If the last message was already from the user, replace it with feedback
    elif messages and messages[-1]["role"] == "user":
        messages[-1]["content"] = feedback
    # Ensure we have the assistant's last response in the history
    elif messages:
        # Add the assistant's last response followed by feedback
        messages.append({"role": "assistant", "content": last_response})
        messages.append({"role": "user", "content": feedback})

    # Ensure we don't exceed token limits by keeping only recent messages if needed
    if len(messages) > 8:  # Arbitrary limit to prevent context length issues
        # Keep system message, last assistant message, and feedback
        system_message = next((m for m in messages if m["role"] == "system"),
                              None)
        new_messages = []
        if system_message:
            new_messages.append(system_message)

        # Add the last few message pairs
        recent_messages = messages[-4:]  # Last 4 messages
        for msg in recent_messages:
            if msg not in new_messages:  # Avoid duplicating system message
                new_messages.append(msg)

        return new_messages

    return messages


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
