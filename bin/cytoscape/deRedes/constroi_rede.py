from typing import List, Dict, Optional
from pathlib import Path
import sys
import os
import gc
import pandas as pd
import py4cytoscape as p4c
import time

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from preprocessing.Processa import ProcessadorDados
from utils.tracer import tracer, trace_step

gc.enable()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
_ESTE_ARQUIVO = Path(__file__).resolve()
RAIZ          = _ESTE_ARQUIVO.parent.parent.parent.parent

OUTPUTS = RAIZ / "results" / "cytoscape" 
OUTPUTS_REDES = OUTPUTS / "redes"
OUTPUTS_TABELAS = OUTPUTS / "tabelas"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

class ConstrutorRede:
    """Classe responsável por orquestrar a construção de redes no Cytoscape."""
    
    def __init__(self, string_cutoff: float = 0.5):
        self.string_cutoff = string_cutoff
        
        # Garantir que diretórios existam
        OUTPUTS_REDES.mkdir(parents=True, exist_ok=True)
        OUTPUTS_TABELAS.mkdir(parents=True, exist_ok=True)
        
    @trace_step("Conectando ao Cytoscape")
    def conectar_cytoscape(self) -> bool:
        """Verifica a conexão com o Cytoscape."""
        try:
            p4c.cytoscape_ping()
            return True
        except Exception:
            return False

    @trace_step("Limpando Nós Desconectados")
    def remover_nos_desconectados(self, network_name: str):
        """Remove nós que não possuem nenhuma aresta."""
        try:
            p4c.set_current_network(network_name)
            nos = p4c.select_nodes_with_degree(0, network=network_name)
            if nos and len(nos) > 0:
                p4c.delete_selected_nodes(network=network_name)
        except Exception as e:
            tracer.console.print(f"    [dim]󰛨 [INFO][/dim] Nenhum nó desconectado removido.")

    @trace_step("Mapeando Estilo Visual")
    def aplica_estilo_visual(self, style_name: str, df_pandas: pd.DataFrame, network_name: str):
        """Aplica o mapeamento visual: cor por log2FC e tamanho por abs_log2FC."""
        try:
            p4c.set_current_network(network_name)
            
            if style_name not in p4c.get_visual_style_names():
                p4c.create_visual_style(style_name)
            
            # 1. Cor do Nó
            v_values = [-7, -5, -2, 0, 2, 5, 7]
            v_colors = ["#503FCF", "#6170C4", "#85A2D2", "#C2DAD9", "#E28383", "#B94A4A", "#AE2020"]
            p4c.set_node_color_default('#D3D3D3', style_name=style_name)
            p4c.set_node_color_mapping('log2FoldChange', v_values, v_colors, mapping_type='c', style_name=style_name)
            
            # 2. Tamanho do Nó
            if 'abs_log2FC' in df_pandas.columns:
                p4c.set_node_size_mapping('abs_log2FC', [df_pandas['abs_log2FC'].min(), df_pandas['abs_log2FC'].max()], [20, 200], mapping_type='c', style_name=style_name)
            
            # 3. Rótulos
            p4c.set_node_label_mapping('display name', style_name=style_name)
            
            # 4. Espessura das Arestas (Busca Dinâmica da Coluna)
            try:
                colunas_edge = p4c.get_table_column_names(table='edge', network=network_name)
                # O STRING app pode nomear como 'experiments', 'stringdb::experiments', etc.
                col_exp = next((c for c in colunas_edge if "experiments" in c.lower()), None)
                
                if col_exp:
                    p4c.set_edge_line_width_mapping(col_exp, [0.0, 0.983], [1, 10], mapping_type='c', style_name=style_name)
                else:
                    tracer.console.print(f"    [dim]󰛨 [INFO][/dim] Coluna de experimentos não encontrada nas arestas.")
            except Exception as e_edge:
                tracer.console.print(f"    [dim]󰛨 [DEBUG][/dim] Falha ao mapear arestas: {e_edge}")

            # Propriedades padrão
            p4c.set_node_font_size_default(14, style_name=style_name)
            p4c.set_edge_color_default('#CCCCCC', style_name=style_name)
            
            p4c.set_visual_style(style_name)
            
        except Exception as e:
            raise e

    @trace_step("Construindo Rede STRING")
    def construir_rede_gene(self, gene_id: str, direction: str, limit: int, df_polars) -> Optional[str]:
        """Constrói a rede para um gene/análise específico."""
        if len(df_polars) > 200:
            df_polars = df_polars.head(200)
            
        gene_string = ProcessadorDados.extrai_genes_string(df_polars)
        if not gene_string:
            return None

        # 1. Query STRING
        string_cmd = (
            f'string protein query query="{gene_string}" '
            f'species="Homo sapiens" '
            f'limit={limit} '
            f'cutoff={self.string_cutoff} '
            f'networkType="full STRING network"'
        )
        
        try:
            # Sub-step visual para a query
            node_q = tracer.push(f"Query STRING (limite {limit}, cutoff {self.string_cutoff})")
            p4c.commands.commands_run(string_cmd)
            time.sleep(8)
            tracer.pop(node_q, f"Query STRING completa", success=True, duration=8.0)
            
            clean_dir = direction.replace('/', '_')
            network_title = f"Rede_{gene_id}_{clean_dir}_lim{limit}"
            p4c.rename_network(network_title)
            
            # 2. Limpeza
            self.remover_nos_desconectados(network_title)
            
            # 3. Mapeamento de Metadados
            node_m = tracer.push("Importando Tabela de Nós e Metadados")
            df_pandas = df_polars.to_pandas()
            for target_col in ['display name', 'query term', 'shared name', 'name']:
                try:
                    col_key = "gene_symbol" if "gene_symbol" in df_pandas.columns else "symbol"
                    p4c.load_table_data(df_pandas, data_key_column=col_key, table_key_column=target_col, network=network_title)
                    break
                except:
                    continue
            tracer.pop(node_m, "Mapeamento de metadados finalizado", success=True, duration=1.0)

            # 4. Estilo
            self.aplica_estilo_visual(f"Estilo_{gene_id}_{clean_dir}", df_pandas, network_title)
            
            # 5. Exportação
            self.exportar_resultados(f"{gene_id}_{clean_dir}_lim{limit}", network_title)
            
            return network_title

        except Exception as e:
            raise e

    @trace_step("Exportando Imagens e Tabelas")
    def exportar_resultados(self, base_name: str, network_title: str):
        analysis_folder = base_name.split("_lim")[0].rsplit("_", 1)[0]
        gene_out_redes = OUTPUTS_REDES / analysis_folder / base_name
        gene_out_tabelas = OUTPUTS_TABELAS / analysis_folder / base_name
        gene_out_redes.mkdir(parents=True, exist_ok=True)
        gene_out_tabelas.mkdir(parents=True, exist_ok=True)

        p4c.set_current_network(network_title)

        node_l = tracer.push("Aplicando Layouts e Gerando Imagens")
        try:
            p4c.layout_network(layout_name='force-directed', network=network_title)
            p4c.export_image(str(gene_out_redes / f"{base_name}_network_force.png"), type='PNG', overwrite_file=True, network=network_title)
            
            p4c.layout_network(layout_name='circular', network=network_title)
            p4c.export_image(str(gene_out_redes / f"{base_name}_network_circular.png"), type='PNG', overwrite_file=True, network=network_title)
        except:
            pass
        tracer.pop(node_l, "Exportação de imagens concluída", success=True, duration=2.0)

        p4c.save_session(str(gene_out_redes / f"{base_name}_session.cys"))
        
        node_table = p4c.get_table_columns(table='node', network=network_title)
        edge_table = p4c.get_table_columns(table='edge', network=network_title)
        node_table.to_csv(gene_out_tabelas / f"{base_name}_string_nodes.csv", index=False)
        edge_table.to_csv(gene_out_tabelas / f"{base_name}_string_edges.csv", index=False)
