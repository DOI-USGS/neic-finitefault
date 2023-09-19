import typer

from .manage import app as manage_app
from .model import app as model_app
from .process import app as process_app

app = typer.Typer(help="CLI for managing the WASP code")
app.add_typer(manage_app, name="manage")
app.add_typer(process_app, name="process-data")
app.add_typer(model_app, name="model")
