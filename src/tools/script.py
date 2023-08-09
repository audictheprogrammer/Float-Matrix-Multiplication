import sys

if __name__ == "__main__":
    filenames = ["../data/benchmark_float_18.txt",
                 "../data/benchmark_float_20.txt",
                 "../data/benchmark_float_22.txt",
                 "../data/benchmark_float_24.txt",
                 "../data/benchmark_float_26.txt",
                 "../data/benchmark_integer_18.txt",
                 "../data/benchmark_integer_20.txt",
                 "../data/benchmark_integer_22.txt",
                 "../data/benchmark_integer_24.txt",
                 "../data/benchmark_integer_26.txt"]
    P = [18, 20, 22, 24, 26];

    text1 = str()
    text2 = str()

    for (i, filename) in enumerate(filenames):
        file = open(filename, "r")
        for line in file:
            temp = line.split();
            if temp[0] == "4096":
                if i < 5:
                    text1 += temp[1] + " " + str(P[i%5]) + "\n"
                else:
                    text2 += temp[1] + " " + str(P[i%5]) + "\n"
        file.close()

    # Float
    file = open("../data/benchmark_recap_float.txt", "w");
    print(text1)
    file.write(text1)
    file.close()
    # Integer
    file = open("../data/benchmark_recap_integer.txt", "w");
    print(text2)
    file.write(text2)
    file.close()
