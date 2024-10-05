# подсчёт GC-состава рида.
def calculate_gc_bounds(seq):
    """
    Выполняет подсчёт Q-score по формуле Unicode-значение символа - 33

    Parameters
    ----------
    seq : str
        последовательность нуклеотидов

    Returns
    -------
    gc_content : float
        Выводит значение содержание "G" и "C" в последовательности
    """
    count_g = seq.lower().count("g")
    count_c = seq.lower().count("c")
    gc_content = ((count_g + count_c) / len(seq)) * 100
    return gc_content


# подсчёт качества рида.
def calculate_quality_threshold(quality):
    """
    Выполняет подсчёт Q-score по формуле Unicode-значение символа - 33

    Parameters
    ----------
    quality : str
        строка из символов и букв, характеризующих
        качество прочтения нуклеотидов

    Returns
    -------
    q_score : float
        Выводит значение среднего Q-score
    """
    seq_length = len(quality)
    summ_q_score = 0
    for letter in quality:
        letter_q_score = ord(letter) - 33
        summ_q_score += letter_q_score

    q_score = summ_q_score / seq_length
    return q_score
