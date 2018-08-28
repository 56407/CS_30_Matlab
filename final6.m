x = [10 8 10 16 15 2 6 7 8 10 12 7 8 9 5 4 2];

nitial_fee = 5000
current_fee =0
current_fee = initial_fee
total_fee = initial_fee
for i in range(3):
current_fee = current_fee + current_fee * increase[i] * .01
total_fee = total_fee + current_fee
print "total fee =",total_fee

Output : total fee = 22974.0