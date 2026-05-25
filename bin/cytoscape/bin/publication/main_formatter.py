import yaml
from pathlib import Path
from table_formatter import AbntTableFormatter
from figure_formatter import AcademicFigureFormatter

# Caminhos Base
_ESTE_ARQUIVO = Path(__file__).resolve()
RAIZ          = _ESTE_ARQUIVO.parent.parent.parent
DOCS_OUT      = RAIZ / "docs" / "manustript"
TABLES_OUT    = DOCS_OUT / "tables"
FIGURES_OUT   = DOCS_OUT / "figures"

def carregar_config():
    config_path = Path(__file__).parent / "config_pub.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def run_publication_pipeline():
    config = carregar_config()
    
    # Garantir diretórios de saída
    TABLES_OUT.mkdir(parents=True, exist_ok=True)
    FIGURES_OUT.mkdir(parents=True, exist_ok=True)
    
    # 1. Inicializar Formatadores
    t_formatter = AbntTableFormatter(
        font_family=config.get('font_family', 'Arial'),
        font_size=config.get('font_size_body', 10)
    )
    
    f_formatter = AcademicFigureFormatter(
        dpi=config.get('dpi', 300),
        font_family=config.get('font_family', 'Arial')
    )
    
    print("--- Iniciando Formatação Acadêmica (ABNT) ---")
    
    # 2. Formatar Tabelas
    print("\n[TABELAS]")
    for t_config in config.get('tables_to_format', []):
        input_p = RAIZ / t_config['input']
        output_p = TABLES_OUT / f"{t_config['output_name']}.xlsx"
        
        if input_p.exists():
            print(f"  > Formatando: {t_config['output_name']}")
            t_formatter.formatar_tabela(input_p, output_p, t_config['title'])
        else:
            print(f"  [AVISO] Arquivo não encontrado: {t_config['input']}")

    # 3. Criar Painéis de Figuras
    print("\n[FIGURAS / PAINÉIS]")
    for f_config in config.get('panels_to_create', []):
        output_p = FIGURES_OUT / f_config['output_name']
        img_paths = [RAIZ / p for p in f_config['images']]
        
        # Filtrar apenas imagens existentes
        existing_imgs = [p for p in img_paths if p.exists()]
        
        if len(existing_imgs) == len(img_paths):
            print(f"  > Criando Painel: {f_config['output_name']}")
            f_formatter.criar_painel(
                existing_imgs, 
                output_p, 
                f_config['labels'], 
                f_config['layout']
            )
            # Salvar legenda em TXT ao lado para facilitar uso com Ollama/Word
            with open(output_p.with_suffix('.txt'), 'w', encoding='utf-8') as f:
                f.write(f_config['caption'])
        else:
            missing = set(img_paths) - set(existing_imgs)
            print(f"  [AVISO] Pulando painel {f_config['output_name']} por falta de imagens: {[p.name for p in missing]}")

    print("\n--- Formatação Concluída! ---")
    print(f"Resultados disponíveis em: {DOCS_OUT}")

if __name__ == "__main__":
    run_publication_pipeline()
