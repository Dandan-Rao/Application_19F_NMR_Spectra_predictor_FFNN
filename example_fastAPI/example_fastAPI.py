from fastapi import FastAPI
from pydantic import BaseModel
from fastapi import HTTPException

# Create an instance of the FastAPI class
app = FastAPI()

# Define a root endpoint
@app.get("/")
def read_root():
    return {"message": "Hello, World!"}

@app.get("/items/{item_id}")
def read_item(item_id: int):
    if item_id == 0:
        raise HTTPException(status_code=404, detail="Item not found")
    return {"item_id": item_id}

class Item(BaseModel):
    name: str
    description: str = None
    price: float
    tax: float = None

@app.post("/items/")
def create_item(item: Item):
    return item

@app.post("/items/", response_model=Item)
def create_item(item: Item):
    return item
