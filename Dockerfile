# Use conda base image
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file and create conda environment
COPY environment.yml .
RUN apt-get update && apt-get install -y build-essential
RUN conda env create -f environment.yml

# Copy application code
COPY src/ ./src/
COPY viewer/ ./viewer/
COPY data/ ./data/
COPY start_backend.sh .
COPY .project-root .

# Create models directory before downloading
RUN mkdir -p aizynthfinder/models/
# Download USPTO models during build
RUN conda run -n deepretro download_public_data aizynthfinder/models/

# Create necessary directories
RUN mkdir -p logs cache_api

# Make start script executable
RUN chmod +x start_backend.sh

# Create non-root user for security
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# Expose port
EXPOSE 5000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD conda run -n deepretro curl -f http://localhost:5000/api/health || exit 1

# Set environment variables
ENV PYTHONPATH=/app
ENV FLASK_APP=src.api
ENV FLASK_ENV=production

# Start the application using conda run (not conda activate)
CMD ["conda", "run", "-n", "deepretro", "python", "src/api.py"] 