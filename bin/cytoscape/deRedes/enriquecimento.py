from pathlib import Path
import py4cytoscape as p4c
import polars as pl
import time
from typing import Optional
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.tracer import trace_step, tracer

class AnalisadorEnriquecimento:
    """Classe responsável por realizar análise de enriquecimento funcional via STRING."""
    
    def __init__(self, base_output_dir: Path):
        self.base_output_dir = base_output_dir
        self.base_output_dir.mkdir(parents=True, exist_ok=True)

    @trace_step("Calculando Enriquecimento STRING")
    def executa_enriquecimento_string(self, network_title: str, gene_name: str) -> Optional[Path]:
        """
        Executa a análise de enriquecimento funcional usando o app STRING no Cytoscape.
        Caso falhe, salva a lista de genes para uso manual.
        """
        analysis_folder = gene_name.split("_lim")[0].rsplit("_", 1)[0]
        gene_out_tabelas = self.base_output_dir / analysis_folder / gene_name
        gene_out_tabelas.mkdir(parents=True, exist_ok=True)
        output_file = gene_out_tabelas / f"{gene_name}_enrichment_raw.csv"
        genes_fallback_file = gene_out_tabelas / f"{gene_name}_genes_para_enriquecimento_manual.txt"

        try:
            p4c.set_current_network(network_title)

            # --- FALLBACK DE SEGURANÇA: Salvar lista de genes da rede ---
            try:
                # Pegar símbolos dos nós da rede atual
                node_data = p4c.tables.get_table_columns(table='node', columns=['display name', 'query term', 'shared name'], network=network_title)
                # Tentar extrair símbolos únicos das colunas mais prováveis
                gene_list = []
                for col in ['display name', 'query term', 'shared name']:
                    if col in node_data.columns:
                        gene_list.extend(node_data[col].dropna().unique().tolist())

                gene_list = sorted(list(set(gene_list)))
                if gene_list:
                    with open(genes_fallback_file, 'w', encoding='utf-8') as f:
                        f.write("\n".join(map(str, gene_list)))
                    tracer.console.print(f"    [dim]󰛨 [INFO][/dim] Lista de {len(gene_list)} genes salva para uso manual.")
            except Exception as e_fallback:
                tracer.console.print(f"    [dim]󰛨 [DEBUG][/dim] Falha ao gerar arquivo de fallback de genes: {e_fallback}")

            # 1. Recuperar Enriquecimento via comando STRING
            p4c.commands.commands_run('string retrieve enrichment')
            time.sleep(6) 

            # 2. Exportar a tabela de enriquecimento dinamicamente
            try:
                tabelas_disponiveis = p4c.get_table_list()
                nome_tabela = next((t for t in tabelas_disponiveis if "enrichment" in t.lower()), "String Enrichment")

                enrichment_data = p4c.get_table_columns(table=nome_tabela, network=network_title)
                if enrichment_data is not None and not enrichment_data.empty:
                    enrichment_data.to_csv(output_file, index=False)
                    return output_file
                else:
                    return None
            except Exception as e_api:
                tracer.console.print(f"    [bold yellow]⚠ [WARNING][/bold yellow] Falha ao capturar tabela de enriquecimento: {e_api}")
                return None

        except Exception as e:
            raise e


    @trace_step("Filtrando Resultados com Polars")
    def filtrar_resultados(self, file_path: Path):
        """Usa Polars para filtrar apenas termos significativos (FDR < 0.05)."""
        if not file_path or not file_path.exists():
            return
            
        try:
            try:
                df = pl.read_csv(file_path)
            except:
                # Fallback para separador alternativo
                df = pl.read_csv(file_path, separator='\t')

            cols = df.columns
            # O STRING pode usar 'FDR' ou 'false discovery rate'
            fdr_col = next((c for c in cols if 'discovery rate' in c.lower() or c.lower() == 'fdr'), None)
            
            if fdr_col:
                # Converter para float se necessário e filtrar
                df_sig = df.filter(pl.col(fdr_col).cast(pl.Float64) < 0.05).sort(fdr_col)
                output_filtered = file_path.with_name(file_path.name.replace("_raw", "_filtered"))
                df_sig.write_csv(output_filtered)
            else:
                tracer.console.print(f"    [bold yellow]⚠ [WARNING][/bold yellow] Coluna FDR não encontrada nas colunas: {cols}")
        except Exception as e:
            tracer.console.print(f"    [bold red]󰅙 Erro ao filtrar enriquecimento:[/bold red] {e}")
