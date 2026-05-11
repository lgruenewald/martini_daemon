from datetime import datetime

date_format = "%Y-%m-%d %H:%M:%S,%f"


def extract_timings_from_log(path: str, ignore_first: bool = False) -> dict[str, float]:
    """
    Given a Martini Daemon log file, extract the timing information of different components.

    This can be used for benchmarking purposes.

    :param path: Path to the log file.
    :param ignore_first: If True, ignore the first frame, skipping the startup cost.
    :return: Dictionary of category to time in seconds.
    """
    categories_starts = {}
    categories_sums = {}
    try:
        with open(path) as f:
            for line in f.readlines():
                date = line[:23]
                content = line[24:]
                parsed_date = datetime.strptime(date, date_format)
                if "start" in content:
                    category = content.replace(" start", "").strip()
                    categories_starts[category] = parsed_date
                if "finished" in content:
                    category = content.replace(" finished", "").strip()
                    prev_start = categories_starts.get(category)
                    if prev_start is None:
                        continue
                    diff = parsed_date - prev_start
                    prev_diff = categories_sums.get(category) or 0.0
                    if ignore_first and categories_sums.get(category) is None:
                        categories_sums[category] = 0.0
                    else:
                        categories_sums[category] = diff.total_seconds() + prev_diff
    except (UnicodeDecodeError, ValueError):
        # not a text file / not a log file with dates
        raise ValueError(f"{path} is not a valid log file.")
    if len(categories_sums) == 0:
        # not a daemon log file we recognize
        raise ValueError(f"{path} is not a valid log file.")
    for category in list(categories_sums.keys()):
        if categories_sums[category] == 0:
            del categories_sums[category]
    return categories_sums
