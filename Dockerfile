# Use AWS Lambda Python 3.12 base image
FROM public.ecr.aws/lambda/python:3.12

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV MPLCONFIGDIR=/tmp
ENV PYTHONPATH=/var/task/external

# Set working directory
WORKDIR /var/task

# Install system dependencies required for scientific computing
RUN microdnf update -y && \
    microdnf install -y \
    gcc \
    gcc-c++ \
    libXrender \
    libXext \
    libSM \
    expat \
    && microdnf clean all

# Copy requirements first for better Docker layer caching
COPY requirements.txt .

# Install Python dependencies with optimizations for Lambda
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --only-binary=all --prefer-binary -r requirements.txt

# Copy application code (including external modules)
COPY . .

# Create necessary directories
RUN mkdir -p /tmp/temp /tmp/artifacts

# Set proper permissions (Lambda runs as lambda user by default)
RUN chmod -R 755 /var/task

# The CMD should specify the handler function
# Format: filename.function_name
CMD ["lambda_handler.lambda_handler"]