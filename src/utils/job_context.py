"""Manages job-specific context using Python's `contextvars`.

This module defines context variables that can be used to store and
access job-specific information throughout an asynchronous call stack or
across different parts of the application handling a specific job or request.

It is particularly useful with structured logging libraries like `structlog`,
where `structlog.contextvars.merge_contextvars` can automatically include
the values of these context variables (e.g., a job ID or a pre-configured
job-specific logger) into log messages.

Attributes:
    logger (contextvars.ContextVar):
A context variable intended to hold job-specific information.
        This could be a job identifier (e.g., a string UUID) or a fully
        configured logger instance that is specific to a job. When set, its
        value can be automatically incorporated into log entries by `structlog`
        if `merge_contextvars` is part of the processor chain.
        The name 'logger' suggests it might hold a logger instance, but it could
        also be used for a job ID string or other contextual data that should
        be logged.
"""
import contextvars

# Create a context variable for the job ID
logger: contextvars.ContextVar = contextvars.ContextVar('logger')
