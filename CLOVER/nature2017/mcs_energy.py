
from decimal import Decimal

lv = 2.26*1e6   # energy of evaporation in J / kg of water
seconds_day = (60*60*24)
storm_acc_rain = 19  # mm / storm on average or kg / m2
wetted_area = 100000*1e6 # storm area in square metres   12 m/s over 24 hours





storm_per_m2 = storm_acc_rain * lv # J / m2 = W s / m2
print('%.2E' % Decimal(storm_per_m2), 'J/m2')


lifetime = storm_per_m2 / seconds_day # W / m2
print('for every m2', lifetime)  # W / m2

watt = lifetime*wetted_area

print('24h storm with 100000km2 wet', '%.2E' % Decimal(watt), 'Watt')

print('Watt hours', '%.2E' % Decimal(watt*24))