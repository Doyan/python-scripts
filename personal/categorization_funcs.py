
def Rubrik(df, contains, category):
    df.category = np.where(df.Rubrik.str.contains(contains), category, df.category)
    return df