"""Convert the Kikwit RData file to CSV.  Usage: python3 convert_kikwit.py in.RData out.csv
Needs: pip install pyreadr"""
import sys
import pyreadr

src, dst = sys.argv[1], sys.argv[2]
df = list(pyreadr.read_r(src).values())[0]
df.to_csv(dst, index=False)
print(df.shape, df[["onset", "death"]].sum().to_dict())
