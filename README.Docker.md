# 19F NMR Spectrum Predictor - Docker Deployment Guide

This guide explains how to build, test, and deploy the 19F NMR Spectrum Predictor application using Docker.

## 🐳 Quick Start

### 1. Build the Docker Image

```bash
# Build the image
docker build -t nmr-predictor:latest .

# Or with a specific tag
docker build -t your-dockerhub-username/nmr-predictor:latest .
```

### 2. Test Locally

```bash
# Run the container
docker run -p 5002:5002 nmr-predictor:latest

# Or use docker-compose
docker-compose up --build
```

### 3. Push to Docker Hub

```bash
# Login to Docker Hub
docker login

# Tag your image
docker tag nmr-predictor:latest your-dockerhub-username/nmr-predictor:latest

# Push to Docker Hub
docker push your-dockerhub-username/nmr-predictor:latest
```

## 🏗️ Building the Image

### Prerequisites

- Docker installed and running
- At least 4GB of available RAM for building
- Stable internet connection for downloading dependencies

### Build Command

```bash
docker build -t nmr-predictor:latest .
```

### Build Options

- **Multi-stage build**: The Dockerfile uses a single-stage build for simplicity
- **Layer caching**: Requirements are installed first for better caching
- **Security**: Runs as non-root user
- **Health checks**: Built-in health check endpoint

## 🚀 Running the Container

### Basic Run

```bash
docker run -p 5002:5002 nmr-predictor:latest
```

### With Environment Variables

```bash
docker run -p 5002:5002 \
  -e SECRET_KEY=your-production-secret-key \
  -e FLASK_ENV=production \
  -e PORT=5002 \
  nmr-predictor:latest
```

### With Volume Mounts

```bash
docker run -p 5002:5002 \
  -v $(pwd)/temp:/app/temp \
  -v $(pwd)/artifacts:/app/artifacts \
  nmr-predictor:latest
```

## 🐙 Using Docker Compose

### Start Services

```bash
docker-compose up --build
```

### Stop Services

```bash
docker-compose down
```

### View Logs

```bash
docker-compose logs -f nmr-predictor
```

## 🔧 Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `FLASK_ENV` | `production` | Flask environment mode |
| `SECRET_KEY` | `your-secret-key-change-in-production` | Flask secret key |
| `PORT` | `5002` | Application port |
| `MPLCONFIGDIR` | `/tmp` | Matplotlib configuration directory |

### Port Configuration

The application runs on port 5002 by default. You can change this by setting the `PORT` environment variable.

## 🏥 Health Checks

The container includes a health check endpoint at `/health` that returns:

```json
{
  "status": "healthy",
  "service": "19F NMR Spectrum Predictor"
}
```

## 📦 Image Details

### Base Image
- **Python 3.11-slim**: Optimized for scientific computing packages

### Size Optimization
- Uses slim base image
- Removes build dependencies after installation
- Excludes unnecessary files via `.dockerignore`

### Security Features
- Runs as non-root user (`appuser`)
- Minimal system packages installed
- Health check endpoint for monitoring

## 🚢 Deployment

### Docker Hub

1. **Build and tag your image**:
   ```bash
   docker build -t your-username/nmr-predictor:latest .
   ```

2. **Push to Docker Hub**:
   ```bash
   docker push your-username/nmr-predictor:latest
   ```

3. **Deploy from other services**:
   ```bash
   docker pull your-username/nmr-predictor:latest
   ```

### Cloud Platforms

#### AWS ECS/Fargate
```yaml
# task-definition.json
{
  "family": "nmr-predictor",
  "containerDefinitions": [
    {
      "name": "nmr-predictor",
      "image": "your-username/nmr-predictor:latest",
      "portMappings": [
        {
          "containerPort": 5002,
          "protocol": "tcp"
        }
      ],
      "environment": [
        {
          "name": "FLASK_ENV",
          "value": "production"
        }
      ]
    }
  ]
}
```

#### Google Cloud Run
```bash
gcloud run deploy nmr-predictor \
  --image your-username/nmr-predictor:latest \
  --platform managed \
  --port 5002 \
  --allow-unauthenticated
```

#### Azure Container Instances
```bash
az container create \
  --resource-group myResourceGroup \
  --name nmr-predictor \
  --image your-username/nmr-predictor:latest \
  --ports 5002 \
  --dns-name-label nmr-predictor
```

## 🔍 Troubleshooting

### Common Issues

1. **Port already in use**:
   ```bash
   # Check what's using port 5002
   lsof -i :5002
   
   # Use a different port
   docker run -p 8080:5002 nmr-predictor:latest
   ```

2. **Permission denied**:
   ```bash
   # Ensure temp directory has correct permissions
   mkdir -p temp && chmod 755 temp
   ```

3. **Memory issues**:
   ```bash
   # Increase Docker memory limit
   # In Docker Desktop: Settings > Resources > Memory
   ```

### Logs

```bash
# View container logs
docker logs <container_id>

# Follow logs in real-time
docker logs -f <container_id>
```

## 📊 Monitoring

### Health Check
```bash
curl http://localhost:5002/health
```

### Application Logs
The application logs all requests and errors. Check container logs for debugging.

## 🔐 Security Considerations

- **Change default secret key** in production
- **Use HTTPS** in production deployments
- **Limit container resources** (CPU, memory)
- **Regular security updates** for base images
- **Network isolation** in production environments

## 📝 License

This Docker configuration is provided as-is for deploying the 19F NMR Spectrum Predictor application.

## 🤝 Support

For issues with Docker deployment:
1. Check the troubleshooting section
2. Review container logs
3. Ensure all prerequisites are met
4. Verify network connectivity and port availability