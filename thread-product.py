import multiprocessing
from itertools import product

alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ12234567890!@#$%^&*?,()-=+[]/;"
num_parts = 4
part_size = len(alphabet) // num_parts

def do_job(first_bits):
    for x in product(first_bits, alphabet, alphabet, alphabet, repeat=22):
        print("".join(x))

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=4)
    results = []
    for i in range(num_parts):
        if i == num_parts - 1:
            first_bit = alphabet[part_size * i :]
        else:
            first_bit = alphabet[part_size * i : part_size * (i+1)]
        results.append(pool.apply_async(do_job(first_bit)))

    pool.close()
    pool.join()