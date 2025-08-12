import click
import requests

API_URL = "http://127.0.0.1:5000/users"  # Base URL of the REST API

@click.group()
def cli():
    """User Management CLI - Interact with the REST API"""
    pass

@click.command()
def get_users():
    """Fetch and display all users"""
    response = requests.get(API_URL)
    if response.status_code == 200:
        users = response.json()
        for user in users:
            click.echo(f"ID: {user['id']}, Name: {user['name']}, Email: {user['email']}")
    else:
        click.echo("Failed to fetch users")

@click.command()
@click.argument("user_id", type=int)
def get_user(user_id):
    """Fetch a single user by ID"""
    response = requests.get(f"{API_URL}/{user_id}")
    if response.status_code == 200:
        user = response.json()
        click.echo(f"ID: {user['id']}, Name: {user['name']}, Email: {user['email']}")
    else:
        click.echo("User not found")

@click.command()
@click.option("--name", prompt="Enter name", help="User's name")
@click.option("--email", prompt="Enter email", help="User's email")
def create_user(name, email):
    """Create a new user"""
    data = {"name": name, "email": email}
    response = requests.post(API_URL, json=data)
    if response.status_code == 201:
        click.echo("User created successfully!")
    else:
        click.echo("Failed to create user")

@click.command()
@click.argument("user_id", type=int)
def delete_user(user_id):
    """Delete a user by ID"""
    response = requests.delete(f"{API_URL}/{user_id}")
    if response.status_code == 200:
        click.echo("User deleted successfully")
    else:
        click.echo("Failed to delete user")

# Add commands to the CLI
cli.add_command(get_users)
cli.add_command(get_user)
cli.add_command(create_user)
cli.add_command(delete_user)

if __name__ == "__main__":
    cli()
