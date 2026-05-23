from typing import List, Dict
from pathlib import Path
import polars as pl
import pandas as pd
import gc

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
_ESTE_ARQUIVO = Path(__file__).resolve()
RAIZ          = _ESTE_ARQUIVO.parent.parent.parent.parent

DADOS = RAIZ / "results"
OUTPUTS = RAIZ / "data" / "processed" / "redes"

LRRK2 = DADOS / "LRRK2_Ranked.xlsx"
SNCA = DADOS / "SNCA_Ranked.xlsx"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def coleta_planilhas(file_path: Path, planilhas: List) -> None:
    dfs_dict = pl.read_excel(
        file_path,
        sheet_name=planilhas,
        engine="calamine"
        )
    for nome, df in dfs_dict.items():
       arquivo_saida = f"{nome}.csv"
       df.write_csv(OUTPUTS / arquivo_saida)
       
       print(f"Planilha '{nome}' processada e salva em: {OUTPUTS / arquivo_saida}")
        
    return 

def concatena_up_down(arquivo_up: Path, arquivo_down: Path) -> pd.DataFrame:
    GENE = arquivo_up.stem.split("_")[0]  # Exemplo: "LRRK2" ou "SNCA"
    print(f"Processando {GENE}...")
    
    lf_up = pl.scan_csv(arquivo_up)
    lf_down = pl.scan_csv(arquivo_down)

    print("Concatenando os DataFrames...")
    q = (
        pl.concat([lf_up, lf_down])
        .with_columns(
            pl.when(pl.col("log2FoldChange") > 0)
            .then(pl.lit("Up"))
            .otherwise(pl.lit("Down"))
            .alias("UP/DOWN")
        )
    )

    q.collect().write_csv(OUTPUTS / f"{GENE}_Up_Down.csv")
        
    return q.collect().to_pandas()

LRRK2_UP = OUTPUTS / "LRRK2_Isogenic_Up.csv"
LRRK2_DOWN = OUTPUTS / "LRRK2_Isogenic_Down.csv"

SNCA_UP = OUTPUTS / "SNCA_Up.csv"
SNCA_DOWN = OUTPUTS / "SNCA_Down.csv"

def main():
    planilhas = ["LRRK2_Isogenic_Up", "LRRK2_Isogenic_Down", "SNCA_Up", "SNCA_Down"]
    coleta_planilhas(LRRK2, planilhas[:2])
    coleta_planilhas(SNCA, planilhas[2:])
    lrrk2_df = concatena_up_down(LRRK2_UP, LRRK2_DOWN)
    snca_df = concatena_up_down(SNCA_UP, SNCA_DOWN)

    print(f"LRRK2 DataFrame:\n{type(lrrk2_df)}\n{lrrk2_df.head()}")
    print(f"SNCA DataFrame:\n{type(snca_df)}\n{snca_df.head()}")

    
    return lrrk2_df, snca_df

if __name__ == "__main__":
    main()
    gc.collect()
