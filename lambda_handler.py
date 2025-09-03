import os
import sys
import json
import logging
import traceback

# Set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Ensure required directories exist
required_dirs = ['/tmp/temp', '/tmp/artifacts']
for dir_name in required_dirs:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=True)

# Global variables for caching
_app = None
_awsgi = None

def get_app():
    """Get Flask app with lazy loading"""
    global _app, _awsgi
    if _app is None:
        try:
            logger.info("Loading Flask app...")
            from app import app as flask_app
            import awsgi
            _app = flask_app
            _awsgi = awsgi
            logger.info("Flask app loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load Flask app: {str(e)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise
    return _app, _awsgi

def lambda_handler(event, context):
    """
    AWS Lambda handler with comprehensive error handling
    """
    try:
        logger.info(f"Lambda function started. Request ID: {context.aws_request_id}")
        logger.info(f"Event keys: {list(event.keys()) if event else 'None'}")
        
        # Ensure event has required fields for awsgi
        if 'queryStringParameters' not in event:
            event['queryStringParameters'] = None
        if 'pathParameters' not in event:
            event['pathParameters'] = None
        if 'multiValueHeaders' not in event:
            event['multiValueHeaders'] = {}
        if 'multiValueQueryStringParameters' not in event:
            event['multiValueQueryStringParameters'] = None
        if 'isBase64Encoded' not in event:
            event['isBase64Encoded'] = False
        
        # Get Flask app (with lazy loading)
        app, awsgi = get_app()
        
        # Call the Flask app through awsgi
        logger.info("Calling awsgi.response...")
        response = awsgi.response(app, event, context, base64_content_types={"image/png"})
        
        logger.info(f"Response generated. Status: {response.get('statusCode', 'Unknown')}")
        return response
        
    except Exception as e:
        logger.error(f"Lambda handler error: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        
        # Return a proper error response with JSON string body
        error_response = {
            'error': 'Internal server error',
            'message': str(e),
            'requestId': context.aws_request_id if context else 'unknown',
            'timestamp': str(context.aws_request_time) if context else 'unknown'
        }
        
        return {
            'statusCode': 500,
            'headers': {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*',
                'Access-Control-Allow-Headers': 'Content-Type',
                'Access-Control-Allow-Methods': 'GET, POST, OPTIONS'
            },
            'body': json.dumps(error_response)
        }