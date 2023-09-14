import typer

from .data import app as data_app

app = typer.Typer(help="CLI for managing the WASP code")
app.add_typer(data_app, name="data")
