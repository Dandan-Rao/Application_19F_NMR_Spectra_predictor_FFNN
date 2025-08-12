#!/usr/bin/env python

# After run in termian: $ ./build_fastapi.py
# Open a web browser and go to the designated URL
# Revise the URL by adding /docs to access the FastAPI documentation.
# Follow the instructions to make a prediction.

from fastapi import FastAPI, HTTPException
from predictor import predictor  # Import predictor function
from pydantic import BaseModel
import uvicorn


app = FastAPI()
class PredictionResponse(BaseModel):
    plot_data: str  # Base64-encoded plot image
    table_data: str  # HTML table
    structure_image_base64: str  # Base64-encoded structure image

class Body(BaseModel):
    txt: str

@app.get("/")
def root():
    return {"message": "Hello World"}

@app.post("/predict")
def predict(body: Body):
    try:
        # Call the predictor function with the SMILES input
        plot_data, table_data, structure_image_base64 = predictor(body.txt)
        # Return the response
        return {
            "plot_data": plot_data,
            "table_data": table_data,
            "structure_image_base64": structure_image_base64,
        }

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Unexpected error: {str(e)}")
    
# Run the FastAPI server with Uvicorn
if __name__ == "__main__":
    uvicorn.run("build_fastapi:app", host="127.0.0.1", port=5000, reload=True)