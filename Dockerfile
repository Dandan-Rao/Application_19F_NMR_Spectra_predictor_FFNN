# syntax=docker/dockerfile:1

FROM eclipse-temurin:23-jre as base

# Set Python version
ARG PYTHON_VERSION=3.12.10

# Environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV MPLCONFIGDIR=/tmp

WORKDIR /app

# Create a non-privileged user
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

# Install required system libraries and Python from source
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    libxrender1 \
    libxext6 \
    libsm6 \
    libexpat1 \
    zlib1g-dev \
    libffi-dev \
    libssl-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    ca-certificates \
    curl \
    && \
    # Download and install Python
    wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz && \
    tar xzf Python-${PYTHON_VERSION}.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure --enable-optimizations && \
    make -j$(nproc) && \
    make altinstall && \
    cd .. && \
    rm -rf Python-${PYTHON_VERSION} Python-${PYTHON_VERSION}.tgz && \
    apt-get remove --purge -y build-essential wget curl && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Make sure python points to the newly installed version
RUN ln -s /usr/local/bin/python3.12 /usr/local/bin/python

# Copy all source code
COPY . .

# Make sure the temp directory is owned by the non-privileged user
RUN mkdir -p /app/temp && chown -R appuser:appuser /app/temp

# Install Python packages
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m ensurepip && \
    python -m pip install --no-cache-dir -r requirements.txt

# Switch to non-privileged user
USER appuser

# Expose the port used by Python app
EXPOSE 5002

# Run the Python application
CMD ["python", "app.py"]
