import typer

from .manage import app as manage_app
from .process import app as process_app

app = typer.Typer(help="CLI for managing the WASP code")
app.add_typer(manage_app, name="manage")
app.add_typer(process_app, name="process-data")
