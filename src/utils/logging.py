"""Logging setup"""
import os
import time
import rootutils

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

# logger = logging.getLogger(__name__)

import logging
from datetime import datetime


def setup_generic_logging():
    logger = logging.getLogger('my_app')
    logger.setLevel(logging.INFO)

    # Create a FileHandler for the daily log file
    today = datetime.now().strftime('%Y-%m-%d')
    daily_log_file = f'app_{today}.log'
    daily_handler = logging.FileHandler(daily_log_file)
    daily_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    daily_handler.setFormatter(daily_formatter)

    # Add the handler to the application logger
    logger.addHandler(daily_handler)
