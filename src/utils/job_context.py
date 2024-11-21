import contextvars

# Create a context variable for the job ID
logger = contextvars.ContextVar('logger')
