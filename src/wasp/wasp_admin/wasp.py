import typer

from .manage import app as manage_app

app = typer.Typer(help="CLI for managing the WASP code")
app.add_typer(manage_app, name="manage")
