#!/usr/bin/env python
from flask_example import Flask, jsonify, request

app = Flask(__name__)

users = [
    {"id": 1, "name": "Alice", "email":"alice@example.com"},
    {"id": 2, "name": "Bob", "email":"bob@example.com"}
]

@app.route('/users', methods = ['GET'])
def get_users():
    return jsonify(users)

@app.route('/users/<int:user_id>', methods = ['GET'])
def get_user(user_id):
    user = next((user for user in users if user['id'] == user_id), None)  # next() returns the next item in an iterator
    if user is None:
        return jsonify({"message": "User not found"}), 404
    return jsonify(user)

@app.route('/users', methods = ['POST'])
def create_user():
    user = request.json   # request.json is a dictionary
    users.append(user)
    return jsonify(user), 201

if __name__ == "__main__":
    app.run(debug = True)