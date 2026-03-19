import pandas as pd
import matplotlib.pyplot as plt


## Таблица интерполированных спектров до оптимизации
df = pd.read_csv("data/spectra.csv")
## Веса, полученные после оптимизации в Mathematica
weights_opt = pd.Series({
    "LED365": 2.19266,
    "LED385wide": 2.17465,
    "LED430": 2.49987,
    "LED480": 1.58427,
    "LED535": 0.977379,
    "LED568wide": 0.851105,
    "LED600": 0.871574,
    "LED635": 0.700088,
    "LED680": 0.650535,
    "LED700": 0.733429,
    "LED740": 0.77227,
    "LED765": 1.14568
})

##Проверка нужных колонок и QE
missing = set(weights_opt.index) - set(df.columns)
if missing:
    raise ValueError(f"Missing columns: {missing}")
if "QE" not in df.columns:
    raise ValueError("Missing QE in the dataframe")

##Колонки после оптимизации, умноженные на QE.
cols = df.columns.intersection(weights_opt.index)
df[cols] = df[cols].mul(weights_opt[cols], axis=1)
df[cols] = df[cols].mul(df["QE"], axis=0)
df.to_csv("data/post_optimized.csv", index=False)




x = df["Wavelength_nm"]

##
spectral_cols = list(cols)

##
plt.figure(figsize=(12, 6))
for col in spectral_cols:
    plt.plot(x, df[col], label=col, linewidth=1)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity (a.u.)")
plt.title("Weighted LED spectra")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# 2) суммарная огибающая = сумма всех спектров
envelope = df[spectral_cols].sum(axis=1)

plt.figure(figsize=(12, 6))
plt.plot(x, envelope, linewidth=2, label="Sum envelope")

plt.ylim(0, 1)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Summed intensity (a.u.)")
plt.title("Summed spectral envelope")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()