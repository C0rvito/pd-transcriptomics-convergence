# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
from matplotlib.pyplot import title
import py4cytoscape as p4c

p4c.cytoscape_ping()       # deve retornar 'You are connected to Cytoscape!'
p4c.cytoscape_version_info()
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<
from rich.console import Console
from rich.panel import Panel



def main():
    console = Console()
    nome_programa = "Cytoscape"
    descricao = f"{p4c.cytoscape_ping()}\n{p4c.cytoscape_version_info()}"
    panel = Panel(f"[blue]{descricao}[/blue]", title=nome_programa)
    console.print(panel)

if __name__ == "__main__":
    main()
