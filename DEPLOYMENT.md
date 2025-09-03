# Deployment Guide

This document provides comprehensive deployment instructions for the 19F NMR Spectrum Predictor application across various cloud platforms.

## 🌐 Live Application

**Production URL**: [https://wqdx9jslij.execute-api.us-east-2.amazonaws.com/prod/](https://wqdx9jslij.execute-api.us-east-2.amazonaws.com/prod/)

## 🚀 AWS Lambda Deployment

### Current Production Setup

The application is deployed on AWS Lambda with API Gateway for serverless, scalable access:

#### Features:
- **Serverless Architecture**: Auto-scaling based on demand
- **Global Access**: Available worldwide via AWS infrastructure
- **Cost-Effective**: Pay only for actual usage
- **High Availability**: Built-in AWS reliability and redundancy

#### Technical Configuration:
- **Runtime**: Python 3.12 on AWS Lambda
- **Memory**: 1024 MB (optimized for ML workloads)
- **Timeout**: 5 minutes (handles cold start + processing)
- **API Gateway**: RESTful API with CORS support
- **Container Image**: Custom Docker image with all dependencies

#### Deployment Steps:

1. **Build Docker Image**:
   ```bash
   docker build -t dandanrao/19fnmr-prediction .
   ```

2. **Push to ECR**:
   ```bash
   # Create ECR repository (if not exists)
   aws ecr create-repository --repository-name 19fnmr-prediction
   
   # Get login token
   aws ecr get-login-password --region us-east-2 | docker login --username AWS --password-stdin <account-id>.dkr.ecr.us-east-2.amazonaws.com
   
   # Tag and push
   docker tag dandanrao/19fnmr-prediction:latest <account-id>.dkr.ecr.us-east-2.amazonaws.com/19fnmr-prediction:latest
   docker push <account-id>.dkr.ecr.us-east-2.amazonaws.com/19fnmr-prediction:latest
   ```

3. **Create Lambda Function**:
   ```bash
   aws lambda create-function \
     --function-name 19fnmr-predictor \
     --package-type Image \
     --code ImageUri=<account-id>.dkr.ecr.us-east-2.amazonaws.com/19fnmr-prediction:latest \
     --role arn:aws:iam::<account-id>:role/lambda-execution-role \
     --timeout 300 \
     --memory-size 1024
   ```

4. **Create API Gateway**:
   ```bash
   # Create REST API
   aws apigateway create-rest-api --name 19fnmr-predictor-api
   
   # Create resource and method
   aws apigateway put-method --rest-api-id <api-id> --resource-id <resource-id> --http-method ANY --authorization-type NONE
   
   # Set up Lambda integration
   aws apigateway put-integration --rest-api-id <api-id> --resource-id <resource-id> --http-method ANY --type AWS_PROXY --integration-http-method POST --uri arn:aws:apigateway:us-east-2:lambda:path/2015-03-31/functions/arn:aws:lambda:us-east-2:<account-id>:function:19fnmr-predictor/invocations
   ```

#### Lambda Configuration:
- **Handler**: `lambda_handler.lambda_handler`
- **Environment Variables**: None required
- **VPC**: Not required (uses default VPC)
- **Dead Letter Queue**: Optional for error handling

## 🐳 Docker Deployment

### AWS ECS/Fargate

```yaml
# task-definition.json
{
  "family": "nmr-predictor",
  "containerDefinitions": [{
    "name": "nmr-predictor",
    "image": "dandanrao/19nmr-predictor_ffnn:latest",
    "portMappings": [{"containerPort": 5002}],
    "healthCheck": {
      "command": ["CMD-SHELL", "curl -f http://localhost:5002/health"],
      "interval": 30,
      "timeout": 5,
      "retries": 3
    },
    "memory": 2048,
    "cpu": 1024
  }]
}
```

### Google Cloud Run

```bash
# Deploy to Cloud Run
gcloud run deploy nmr-predictor \
  --image dandanrao/19nmr-predictor_ffnn:latest \
  --platform managed \
  --port 5002 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300
```

### Azure Container Instances

```bash
# Deploy to Azure
az container create \
  --resource-group myResourceGroup \
  --name nmr-predictor \
  --image dandanrao/19nmr-predictor_ffnn:latest \
  --ports 5002 \
  --dns-name-label nmr-predictor \
  --memory 2 \
  --cpu 2
```

## 🔧 Environment Configuration

### Required Environment Variables

For container deployments:
```bash
FLASK_ENV=production
SECRET_KEY=your-secret-key-here
PORT=5002
```

### Health Check Endpoint

All deployments include a health check endpoint:
```bash
curl http://your-domain/health
```

Expected response:
```json
{
  "status": "healthy",
  "service": "19F NMR Spectrum Predictor"
}
```

## 📊 Performance Optimization

### Lambda Optimization
- **Memory**: 1024 MB (based on actual usage of ~896 MB)
- **Timeout**: 5 minutes (handles 30-second cold start)
- **Provisioned Concurrency**: Consider for production to reduce cold starts

### Container Optimization
- **Memory**: 2 GB minimum for ML workloads
- **CPU**: 2 cores recommended
- **Health Checks**: Configured for reliability

## 🔒 Security Considerations

### Lambda Security
- **IAM Role**: Minimal permissions for Lambda execution
- **VPC**: Not required (uses default VPC)
- **API Gateway**: CORS enabled for web access

### Container Security
- **Non-root User**: Container runs as non-privileged user
- **Image Scanning**: Regular security scans
- **Network Policies**: Restrictive network access

## 🐛 Troubleshooting

### Common Issues

**Lambda Cold Start**:
- First request takes ~30 seconds
- Subsequent requests are fast
- Consider provisioned concurrency for production

**Memory Issues**:
- Monitor CloudWatch metrics
- Increase memory if needed
- Check for memory leaks

**API Gateway Errors**:
- Check Lambda function logs
- Verify CORS configuration
- Test Lambda function directly

### Monitoring

**CloudWatch Metrics**:
- Duration
- Memory usage
- Error rate
- Invocation count

**Logs**:
```bash
# View Lambda logs
aws logs describe-log-groups --log-group-name-prefix /aws/lambda/19fnmr-predictor

# View recent logs
aws logs filter-log-events --log-group-name /aws/lambda/19fnmr-predictor --start-time $(date -d '1 hour ago' +%s)000
```

## 💰 Cost Optimization

### Lambda Costs
- **Memory**: 1024 MB (optimized based on usage)
- **Duration**: ~10 seconds average
- **Requests**: Pay per invocation

### Container Costs
- **Always-on**: Consider for high-traffic applications
- **Auto-scaling**: Configure based on demand
- **Resource limits**: Set appropriate CPU/memory limits

## 🔄 CI/CD Pipeline

### GitHub Actions Example

```yaml
name: Deploy to AWS Lambda

on:
  push:
    branches: [main]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      
      - name: Configure AWS credentials
        uses: aws-actions/configure-aws-credentials@v1
        with:
          aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}
          aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
          aws-region: us-east-2
      
      - name: Build and push Docker image
        run: |
          docker build -t $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG .
          docker push $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG
      
      - name: Update Lambda function
        run: |
          aws lambda update-function-code \
            --function-name 19fnmr-predictor \
            --image-uri $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG
```

## 📞 Support

For deployment issues:
- Check CloudWatch logs
- Review Lambda function metrics
- Test health check endpoint
- Verify API Gateway configuration

---

**Note**: This deployment guide covers the most common scenarios. Adjust configurations based on your specific requirements and traffic patterns.
