#!/usr/bin/env python
import os
import logging
from flask import Flask, request, render_template
from predictor import predictor

# import awsgi

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Production configuration
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'your-secret-key-change-in-production')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        smiles = request.form.get("smiles")

        if not smiles or not isinstance(smiles, str):
            logger.warning(f"Invalid SMILES input: {smiles}")
            return render_template(
                "index.html", error_message="Please provide a valid SMILES string."
            )

        try:
            logger.info(f"Processing SMILES: {smiles}")
            
            # Check if required directories exist (use /tmp for writable directories)
            import os
            required_dirs = ['/tmp/temp', '/tmp/artifacts']
            for dir_name in required_dirs:
                if not os.path.exists(dir_name):
                    os.makedirs(dir_name, exist_ok=True)
                    logger.info(f"Created directory: {dir_name}")
            
            # Call the predictor function with the SMILES input
            plot_data, table_data, structure_image_base64 = predictor(smiles)
            
            if plot_data is None or table_data is None:
                raise ValueError("Failed to generate prediction results")
                
            logger.info(f"Successfully processed SMILES: {smiles}")
            return render_template(
                "index.html",
                plot_data=plot_data,
                table_data=table_data,
                structure_image_base64=structure_image_base64,
            )
        except ValueError as e:
            logger.error(f"ValueError for SMILES {smiles}: {str(e)}")
            return render_template("index.html", error_message=str(e))
        except ImportError as e:
            logger.error(f"ImportError for SMILES {smiles}: {str(e)}")
            return render_template("index.html", error_message=f"Missing dependency: {str(e)}")
        except FileNotFoundError as e:
            logger.error(f"FileNotFoundError for SMILES {smiles}: {str(e)}")
            return render_template("index.html", error_message=f"Required file not found: {str(e)}")
        except MemoryError as e:
            logger.error(f"MemoryError for SMILES {smiles}: {str(e)}")
            return render_template("index.html", error_message="Insufficient memory to process the request. Please try a simpler compound.")
        except Exception as e:
            logger.error(f"Unexpected error for SMILES {smiles}: {str(e)}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return render_template(
                "index.html", 
                error_message=f"An unexpected error occurred: {str(e)}. Please try again or contact support."
            )

    return render_template("index.html")

@app.route("/health")
def health_check():
    """Health check endpoint for Docker health checks"""
    return {"status": "healthy", "service": "19F NMR Spectrum Predictor"}

@app.errorhandler(404)
def not_found(error):
    return render_template("index.html", error_message="Page not found"), 404

@app.errorhandler(500)
def internal_error(error):
    logger.error(f"Internal server error: {error}")
    import traceback
    logger.error(f"Traceback: {traceback.format_exc()}")
    return render_template("index.html", error_message="Internal server error. Please try again."), 500




# def lambda_handler(event, context):
#     """AWS Lambda handler using aws-wsgi"""
#     return awsgi.response(app, event, context, base64_content_types={"image/png"})

if __name__ == "__main__":
    # Production configuration
    port = int(os.environ.get('PORT', 5002))
    debug = os.environ.get('FLASK_ENV') == 'development'
    
    logger.info(f"Starting 19F NMR Spectrum Predictor on port {port}")
    app.run(
        host='0.0.0.0',
        port=port,
        debug=debug
    )