
def time_mean(ds, month=True, day=True, year=True, t=True):
    """
    Tries to do the time mean. If no time is detected, return the array unchanged
    """
    if month:
        try:
            ds = ds.mean('month')
        except ValueError:
            pass
    if day:
        try:
            ds = ds.mean('day')
        except ValueError:
            pass
    if year:
        try:
            ds = ds.mean('year')
        except ValueError:
            pass
    if t:
        try:
            ds = ds.mean('t')
        except ValueError:
            pass
    return ds
