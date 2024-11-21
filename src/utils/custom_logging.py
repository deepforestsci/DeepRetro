"""Logging setup"""
import os
import time
import rootutils
from datetime import datetime
from structlog.processors import ProcessorFormatter
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
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.FileHandler(
                f"app_{datetime.now().strftime('%Y-%m-%d')}.log"),
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
            structlog.dev.ConsoleRenderer(),
        ],
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.make_filtering_bound_logger(logging.INFO),
        cache_logger_on_first_use=True,
    )


class JobSpecificFileHandler(logging.FileHandler):

    def __init__(self, job_id, *args, **kwargs):
        filename = f'job_{job_id}.log'
        super().__init__(filename, *args, **kwargs)


def add_job_specific_handler(logger, job_id):
    handler = JobSpecificFileHandler(job_id)
    formatter = ProcessorFormatter(processor=structlog.dev.ConsoleRenderer(), )
    handler.setFormatter(formatter)
    logger._logger.addHandler(handler)
    return handler
