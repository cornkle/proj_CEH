import numpy as np
import matplotlib.pyplot as plt
import random
from skimage.draw import disk
import ipdb
from land_wavelet import wclass


def draw_circle(x, y, radius):

    xloc1 = np.arange(x - radius, x + radius + 1)
    yloc1 = np.arange(y - radius, y + radius + 1)
    xloc, yloc = np.meshgrid(xloc1, yloc1)
    distloc = ( ((xloc - x) * (xloc - x)) + ((yloc - y) * (yloc - y)) )**.5

    indloc = (distloc <= radius).nonzero()
    ycirc = indloc[0] - radius + y
    xcirc = indloc[1] - radius + x

    return (ycirc, xcirc)


def generate_features(image_size, num_features, feature_sizes, clustered=False, fixed_centers=None, remove_largest=False, isolate_features=None, additional_isolated=None, additional_feature_size=None):
    img = np.zeros((image_size, image_size))

    if fixed_centers is not None:
        centers = fixed_centers
    elif clustered:
        clusters = [(int(image_size / 4), int(image_size / 4)),
                    (int(3 * image_size / 4), int(3 * image_size / 4))]
        centers = []
        for _ in range(num_features):
            cluster = random.choice(clusters)
            x = np.clip(int(np.random.normal(cluster[0], image_size / 8)), 0, image_size - 1)
            y = np.clip(int(np.random.normal(cluster[1], image_size / 8)), 0, image_size - 1)
            centers.append((x, y))
    else:
        centers = [tuple(np.random.randint(0, image_size, 2)) for _ in range(num_features)]

    if remove_largest:
        largest_index = np.argmax(feature_sizes)
        feature_sizes = np.delete(feature_sizes, largest_index)
        centers.pop(largest_index)

    if isolate_features is not None:
        for index in isolate_features:
            centers[index] = tuple(np.random.randint(0, image_size, 2))

    if additional_isolated is not None and additional_feature_size is not None:
        centers.extend(additional_isolated)
        feature_sizes = np.append(feature_sizes, [additional_feature_size] * len(additional_isolated))

    for (x, y), size in zip(centers, feature_sizes):
        rr, cc = disk((x, y), size, shape=img.shape)
        img[rr, cc] = 1

    return img, centers


def generate_equidistant_features(image_size, num_features, feature_size):
    img = np.zeros((image_size, image_size))
    spacing = int(np.sqrt((image_size ** 2) / num_features))
    centers = [(x, y) for x in range(spacing // 2, image_size, spacing) for y in range(spacing // 2, image_size, spacing)]
    centers = centers[:num_features]

    for (x, y) in centers:
        rr, cc = disk((x, y), feature_size, shape=img.shape)
        img[rr, cc] = 1

    return img, centers


def plot_images(images, titles, size):
    fig, axs = plt.subplots(2, 4, figsize=(20, 10))
    tick_positions = np.linspace(-size//2, size//2, 5, dtype=int)

    for ax, img, title in zip(axs.ravel(), images, titles):
        cax = ax.imshow(img, cmap='gray_r', extent=[-size//2, size//2, -size//2, size//2])
        ax.set_title(title)
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xlabel('Distance from center (km)')
        ax.set_ylabel('Distance from center (km)')
        ax.grid(True)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(cax, cax=cbar_ax, orientation='vertical', label='Feature Presence')
    plt.suptitle(f'Domain Size: {size}x{size} km')
    plt.tight_layout(rect=[0, 0.03, 0.9, 0.95])
    plt.show()
    fig.savefig('/home/ck/DIR/cornkle/figs/COCOON/organisation/test_cases.png')


def main():
    image_size = 100
    num_features = 10
    feature_size = 5
    random_seed = 66

    random.seed(random_seed)
    np.random.seed(random_seed)

    # Case 1: Equidistant features, same scale
    equidistant_img, _ = generate_equidistant_features(image_size, num_features, feature_size)

    # Case 2: Random locations, same scale
    random_img, random_centers = generate_features(image_size, num_features, [feature_size] * num_features, clustered=False)

    # Case 3: 2 small isolated features
    small_isolated_centers = [(int(image_size * 0.2), int(image_size * 0.2)), (int(image_size * 0.8), int(image_size * 0.8))]
    small_isolated_img, _ = generate_features(image_size, 2, [feature_size] * 2, fixed_centers=small_isolated_centers)

    # Case 4: 2 large isolated features (same centers as Case 3)
    large_isolated_img, _ = generate_features(image_size, 2, [feature_size * 3] * 2, fixed_centers=small_isolated_centers)
    circ_ind1 = draw_circle(50, 50, 5)
    large_isolated_img[circ_ind1] = 1

    # Case 5: Clustered, same scale
    clustered_img, clustered_centers = generate_features(image_size, num_features, [feature_size] * num_features, clustered=True)

    # Case 6: Clustered, 3 different scales (including small scale used in 5)
    different_sizes = [feature_size, feature_size * 2, feature_size * 3]
    feature_sizes_different = np.random.choice(different_sizes, num_features, replace=True)
    clustered_diff_size_img, _ = generate_features(image_size, num_features, feature_sizes_different, clustered=True, fixed_centers=clustered_centers)

    # Case 7: Same as 6, but largest feature removed
    clustered_diff_size_removed_img, _ = generate_features(image_size, num_features, feature_sizes_different, clustered=True, fixed_centers=clustered_centers, remove_largest=True)

    # Case 8: Same as 6, but 2 small features displaced to be isolated instead of clustered
    #small_indices = [i for i, size in enumerate(feature_sizes_different) if size == feature_size][:2]

    #clustered_isolated_img, _ = generate_features(image_size, num_features, feature_sizes_different, clustered=True, fixed_centers=clustered_centers, remove_largest=True)
    circ_ind1 = draw_circle(10,50,5)
    circ_ind2 = draw_circle(90,50,5)

    clustered_isolated_img = clustered_diff_size_removed_img.copy()

    #clustered_isolated_img[circ_ind1] = 1
    clustered_isolated_img[circ_ind2] = 1


    images = [equidistant_img, random_img, small_isolated_img, large_isolated_img, clustered_img, clustered_diff_size_img, clustered_diff_size_removed_img, clustered_isolated_img]

    titles = [
        "Equidistant",
        "Random Locations",
        "2 Small Features",
        "2 Large Features",
        "Clustered, small",
        "Clustered, 3 Scales",
        "Clustered, 3 Scales, Largest Removed",
        "Clustered, 3 Scales, 2 Isolated"
    ]

    strc = {title: img for title, img in zip(titles, images)}
    plot_images(images, titles, image_size)

    return strc



def run_wav(strc):

    wObj = wclass.landwav('organisation')

    p_strc = []
    powers = []
    pure = []


    for k in strc.keys():
        img = strc[k]
        wObj.read_img(img, np.arange(img.shape[1]), np.arange(img.shape[0]))
        coeffs, power, scales, period = wObj.applyWavelet(normed='none', le_thresh=0, fill=0.01)


        frac_power = (np.sum(power,axis=(1,2))/np.sum(power))*100


        p_strc.append(frac_power)
        powers.append(np.sum(power, axis=(1,2)))
        pure.append(np.sum(power, axis=0)) #/np.sum(power)

    print(scales)
    titles = list(strc.keys())

    #ipdb.set_trace()
    f = plt.figure(figsize=(9,5), dpi=300)
    ax = f.add_subplot(111)
    for ids, ps in enumerate(p_strc):
        ax.plot(scales, ps, label=titles[ids])
    ax.legend()
    ax.set_xlabel('Scales (km)')
    ax.set_ylabel('Scale contribution to total power (%)')
    plt.show()
    plt.tight_layout()
    f.savefig('/home/ck/DIR/cornkle/figs/COCOON/organisation/scale_power_contribution.png')


    f = plt.figure(figsize=(9,5), dpi=300)
    ax = f.add_subplot(111)

    for ids, ps in enumerate(powers):
        cdf = np.cumsum(np.array(ps))

        # Normalize the cumulative sum to the range [0, 1]
        cdf_normalized = cdf / cdf[-1]  # Dividing by the last value to get the cumulative percentage

        ax.plot(scales, cdf_normalized, label=titles[ids])

    ax.set_xlabel('Scales (km)')
    ax.set_ylabel('Scale CDF')
    ax.legend()
    plt.show()
    plt.tight_layout()
    f.savefig('/home/ck/DIR/cornkle/figs/COCOON/organisation/scale_CDF.png')

    plot_images(pure, titles, 100)