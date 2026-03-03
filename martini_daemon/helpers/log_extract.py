from datetime import datetime

date_format = "%Y-%m-%d %H:%M:%S,%f"

def extract(path, ignore_first=False):
    categories_starts = {}
    categories_sums = {}
    try:
        with open(path, "r") as f:
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
                        categories_sums[category] = 0.
                    else:
                        categories_sums[category] = diff.total_seconds() + prev_diff
    except (UnicodeDecodeError, ValueError):
        # not a text file / not a log file with dates
        return False
    if len(categories_sums) == 0:
        # not a daemon log file we recognize
        return False
    return categories_sums

