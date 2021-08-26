import re

LIST_ENG = """
hokkaido
aomori
iwate
miyagi
akita
yamagata
hukushima
ibaraki
tochigi
gunma
saitama
chiba
tokyo
kanagawa
nigata
toyama
ishikawa
hukui
yamanashi
nagano
gihu
shizuoka
aichi
mie
shiga
kyoto
osaka
hyogo
nara
wakayama
tottori
shimane
okayama
hiroshima
yamaguchi
tokushima
kagawa
ehime
kochi
hukuoka
saga
nagasaki
kumamoto
oita
miyazaki
kagoshima
okinawa
"""

PREFECTURE_ENG = re.findall(r'(\w+)\n', LIST_ENG)


