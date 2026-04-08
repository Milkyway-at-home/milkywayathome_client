import matplotlib.pyplot as plt

def make_histogram(input_files, output_file="both.png", bins=200, limit=50000, colors=None, labels=None):
    if isinstance(input_files, str):
        input_files = [input_files]
    if colors is None:
        colors = ["steelblue", "tomato", "green", "orange", "purple"]
    if labels is None:
        labels = input_files

    plt.figure()
    for i, input_file in enumerate(input_files):
        data = []
        with open(input_file, "r") as f:
            for _, line in zip(range(limit), f):
                line = line.strip()
                if line:
                    data.append(float(line))
        plt.hist(data, bins=bins, color=colors[i % len(colors)], label=labels[i],
                 alpha=0.4, histtype="stepfilled", edgecolor=colors[i % len(colors)])

    plt.xlabel("Value")
    plt.yscale("log")
    plt.ylabel("Frequency")
    plt.title("Probability")
    plt.legend()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    make_histogram(["halo_prob.txt", "halo_prob.txt"])
