import pandas as pd 
import polars as pl  
import time

# Pandas
start = time.time()
df_pandas = pd.read_csv('dados_gerados.csv')
pandas_time = time.time() - start

# Polars
start = time.time()
df_polars = pl.read_csv('dados_gerados.csv')
polars_time = time.time() - start

print(f"Pandas read time: {pandas_time:.2f} seconds")
print(f"Polars read time: {polars_time:.2f} seconds")
print(f"Polars is {pandas_time/polars_time:.1f}x faster")