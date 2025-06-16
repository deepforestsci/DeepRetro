"""Custom logging setup for the DeepRetro application.

This module configures both standard Python logging and `structlog` to provide
structured, context-aware logging. It establishes:
- Date-based log directories (e.g., `logs/YYYY-MM-DD/`).
- A general application log file (e.g., `app_YYYY-MM-DD.log`) within the
  date-specific directory.
- A stream handler for console output, ensuring logs are visible during
  interactive sessions or in container logs.
- `structlog` processors for essential features like timestamping, log level
  inclusion, call site information (filename and line number), and automatic
  formatting of exception information.
- A mechanism to add job-specific log file handlers. This allows logs related
  to a particular task, request, or "job" (identified by a `job_id`) to be
  written to a dedicated file (e.g., `logs/YYYY-MM-DD/job_<job_id>.log`).

Usage:
    Call `setup_logging()` once at application startup. Afterwards, obtain logger
    instances using `structlog.get_logger()`. To add a job-specific handler,
    use `add_job_specific_handler(logger, job_id)`.
"""
import os
import time
import rootutils
from datetime import datetime
from structlog.stdlib import ProcessorFormatter
import structlog

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)
import logging

# Set up logging
# Create logs directory if it doesn't exist
if not os.path.exists(f'{root_dir}/logs'):
    os.makedirs(f'{root_dir}/logs')

# Create date-based logs directory if it doesn't exist
date_dir = f'{root_dir}/logs/{time.strftime("%Y-%m-%d")}'
if not os.path.exists(date_dir):
    os.makedirs(date_dir)


def setup_logging():
    """Initializes and configures global logging settings.

    This function sets up the basic Python logging handlers (a file handler
    for the main application log and a stream handler for console output).
    It also configures `structlog` with a series of processors to enable
    structured logging with features like timestamps, log levels, and call
    site parameters. This function should typically be called once at the
    beginning of the application's lifecycle.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.FileHandler(
                f"{date_dir}/app_{datetime.now().strftime('%Y-%m-%d')}.log"),
            logging.StreamHandler()
        ])

    structlog.configure(
        processors=[
            structlog.contextvars.merge_contextvars,
            structlog.processors.TimeStamper(fmt='%Y-%m-%d %H:%M:%S'),
            structlog.processors.add_log_level,
            structlog.processors.CallsiteParameterAdder([
                structlog.processors.CallsiteParameter.FILENAME,
                structlog.processors.CallsiteParameter.LINENO
            ]),
            structlog.processors.format_exc_info,
            structlog.dev.ConsoleRenderer(colors=False),
        ],
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.make_filtering_bound_logger(logging.INFO),
        cache_logger_on_first_use=True,
    )


class JobSpecificFileHandler(logging.FileHandler):
    """A logging handler that directs messages to a job-specific log file.

    This handler extends `logging.FileHandler`. When instantiated, it creates
    or appends to a log file named `job_<job_id>.log` located within the
    application's date-specific log directory (e.g., `logs/YYYY-MM-DD/`).

    Args:
        job_id (str): The unique identifier for the job. This ID is used
            to name the log file.
        *args: Additional positional arguments to pass to the constructor
            of `logging.FileHandler`.
        **kwargs: Additional keyword arguments to pass to the constructor
            of `logging.FileHandler`.
    """

    def __init__(self, job_id, *args, **kwargs):
        filename = f'{date_dir}/job_{job_id}.log'
        super().__init__(filename, *args, **kwargs)


def add_job_specific_handler(logger, job_id):
    """Adds a job-specific file handler to a structlog logger instance.

    This function creates an instance of `JobSpecificFileHandler` for the given
    `job_id` and adds it to the underlying standard Python logger of the provided
    `structlog` logger. This facilitates routing logs for a specific job to its
    own dedicated file.

    Args:
        logger (structlog.stdlib.BoundLogger): The `structlog` logger instance
            to which the job-specific handler will be added. It's assumed this
            logger uses the standard library's logging system as its backend.
        job_id (str): The unique identifier for the job, used to create the
            `JobSpecificFileHandler` and name its log file.

    Returns:
        JobSpecificFileHandler: The created and added job-specific handler instance.
            This can be useful if the handler needs to be removed or manipulated later.
    """
    handler = JobSpecificFileHandler(job_id)
    formatter = ProcessorFormatter(
        processor=structlog.dev.ConsoleRenderer(colors=False), )
    handler.setFormatter(formatter)
    logger._logger.addHandler(handler)
    return handler
