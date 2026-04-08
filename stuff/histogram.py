import matplotlib.pyplot as plt

def make_histogram(input_file, output_file="histogram.png", bins=200, limit=100000):
    # Read up to `limit` floats from the file
    data = []
    with open(input_file, "r") as f:
        for _, line in zip(range(limit), f):
            line = line.strip()
            if line:
                data.append(float(line))

    # Plot histogram
    plt.figure()
    plt.hist(data, bins=bins)
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.title(f"Probability")
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    make_histogram("data.txt")
