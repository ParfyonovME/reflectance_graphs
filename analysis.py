import pandas as pd

df = pd.read_csv("data/post_optimized.csv")

# спектральные колонки (исключаем wavelength и QE)
exclude = {"Wavelength_nm", "QE"}
spectral_cols = [c for c in df.columns if c not in exclude]

# суммарный спектр
S = df[spectral_cols].sum(axis=1)

mean = S.mean()
std = S.std()
cv = std / mean
p2p_rel = (S.max() - S.min()) / mean
r_max_min = S.max() / S.min()

## Вывод метрик
print(f"CV: {cv*100:.2f}%")
print(f"P2P_rel: {p2p_rel*100:.2f}%")
print(f"R_max/min: {r_max_min:.3f}")

##Значения слегка отличаются от статейных из-за неотсчения излишнего диапазона справа