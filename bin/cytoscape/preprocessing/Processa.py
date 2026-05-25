from typing import List, Dict
from pathlib import Path
import polars as pl
import gc
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.tracer import trace_step, tracer

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
_ESTE_ARQUIVO = Path(__file__).resolve()
RAIZ          = _ESTE_ARQUIVO.parent.parent.parent.parent
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

class ProcessadorDados:
    """Classe responsável pelo processamento de dados usando Polars."""
    
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    @trace_step("Extraindo Planilhas Excel")
    def coleta_planilhas(self, file_path: Path, planilhas: List[str]) -> Dict[str, Path]:
        """Lê planilhas de um arquivo Excel e salva como CSVs separados."""
        dfs_dict = pl.read_excel(
            file_path,
            sheet_name=planilhas,
            engine="calamine"
        )
        
        caminhos_saida = {}
        for nome, df in dfs_dict.items():
            arquivo_saida = self.output_dir / f"{nome}.csv"
            df.write_csv(arquivo_saida)
            caminhos_saida[nome] = arquivo_saida
        
        del dfs_dict
        gc.collect()
        return caminhos_saida

    def preparar_direcao(self, df_base: pl.LazyFrame, direcao: str) -> pl.DataFrame:
        """Aplica filtros e colunas extras conforme a direção desejada."""
        q = df_base.with_columns(
            pl.when(pl.col("log2FoldChange") > 0)
            .then(pl.lit("Up"))
            .otherwise(pl.lit("Down"))
            .alias("UP/DOWN"),
            pl.col("log2FoldChange").abs().alias("abs_log2FC")
        ).filter(pl.col("padj") < 0.05)

        if direcao == "UP":
            q = q.filter(pl.col("log2FoldChange") > 0)
        elif direcao == "DOWN":
            q = q.filter(pl.col("log2FoldChange") < 0)
        
        return q.sort("padj").collect()

    @trace_step("Processando Análise")
    def processar_analise(self, analysis_id: str, excel_path: Path, up_sheet: str, down_sheet: str, directions: List[str]) -> Dict[str, pl.DataFrame]:
        """Processa uma análise específica e retorna DataFrames para cada direção."""
        planilhas = [up_sheet, down_sheet]
        caminhos = self.coleta_planilhas(excel_path, planilhas)
        
        lf_up = pl.scan_csv(caminhos[up_sheet])
        lf_down = pl.scan_csv(caminhos[down_sheet])
        lf_combined = pl.concat([lf_up, lf_down])

        resultados = {}
        for dir_name in directions:
            if dir_name == "UP":
                df = self.preparar_direcao(lf_up, "UP")
            elif dir_name == "DOWN":
                df = self.preparar_direcao(lf_down, "DOWN")
            else: # UP/DOWN
                df = self.preparar_direcao(lf_combined, "UP/DOWN")
            
            # Salvar arquivo intermediário
            nome_arquivo = f"{analysis_id}_{dir_name.replace('/', '_')}.csv"
            df.write_csv(self.output_dir / nome_arquivo)
            resultados[dir_name] = df
            
        gc.collect()
        return resultados

    @staticmethod
    def extrai_genes_string(df: pl.DataFrame) -> str:
        """Extrai os símbolos dos genes como uma string separada por vírgulas."""
        col_symbol = "gene_symbol" if "gene_symbol" in df.columns else "symbol"
        if col_symbol not in df.columns:
             return ""
        genes = (
            df.select(col_symbol)
            .drop_nulls()
            .unique()
            .to_series()
            .to_list()
        )
        return ",".join(map(str, genes))
