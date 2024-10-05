DIFF_DNA_RNA = {
    "T": "U",
    "t": "u",
    "U": "T",
    "u": "t"
}
DNA_TO_RNA = {
    "A": "U",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "t": "a",
    "g": "c",
    "c": "g",
}
RNA_TO_RNA = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
}
DNA_TO_DNA = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
}
RNA_TO_DNA = {
    "A": "T",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "u": "a",
    "g": "c",
    "c": "g",
}


# transcribe — вернуть транскрибированную последовательность
def transcribe(arg):
    """
    Функция принимает на вход последовательность из ДНК или РНК и возвращает
    транскрибированную последовательность

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность: ДНК или РНК

    Returns
    -------
    result : str
        Транскрибированная последовательность нуклеотидов РНК или ДНК
    """
    result = ""

    for letter in arg:
        if letter in DIFF_DNA_RNA:
            result += DIFF_DNA_RNA[letter]
        else:
            result += letter
    return result


# reverse — вернуть развёрнутую последовательность
def reverse(arg):
    """
    Функция принимает на вход последовательность из ДНК или РНК и возвращает
    развёрнутую последовательность

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность: ДНК или РНК

    Returns
    -------
    result : str
        Развёрнутая последовательность нуклеотидов
    """
    return arg[::-1]


# complement — вернуть комплементарную последовательность
def complement(arg):
    """
    Функция принимает на вход последовательность из ДНК или РНК и возвращает
    комплементарную ей последовательность

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность: ДНК или РНК

    Returns
    -------
    result : str
        Комплементарная последовательность нуклеотидов
    """
    result = ""
    for letter in arg:
        if letter in DNA_TO_DNA:
            result += DNA_TO_DNA[letter]
        else:
            result += RNA_TO_RNA[letter]
    return result


# reverse_complement — вернуть обратную комплементарную последовательность
def reverse_complement(arg):
    """
    Функция принимает на вход последовательность из ДНК или РНК и возвращает
    комплементарную и развёрнутую последовательность

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность: ДНК или РНК

    Returns
    -------
    result : str
        Развёрнутая комплементарная последовательность нуклеотидов
    """
    arg = complement(arg)
    return reverse(arg)


# is_palindrome — является ли поданная последовательность биопалиндромом
def is_palindrome(arg) -> bool:
    """
    Проверяет, соответствует ли поданная последовательность
    комплементарной развёрнутой последовательности

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность

    Returns
    -------
    bool
        True если является биопалиндромом, False если не является
    """
    reverse = ""
    for letter in arg.lower():
        if letter in DNA_TO_DNA:
            reverse += DNA_TO_DNA[letter]
        else:
            reverse += RNA_TO_RNA[letter]

    return arg.lower() == reverse[::-1].lower()


# annealing_temperature — подсчёт температуры отжига праймера
def annealing_temperature(arg):
    """
    Считает температуру отжига последовательности по формуле
    Tm (°C ) = 2 х (A+T) + 4 х (G+C)

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность

    Returns
    -------
    result : float
        Температура отжига последовательности
    """
    arg = arg.lower()
    count_a = arg.count("a")
    count_t = arg.count("t")
    count_g = arg.count("g")
    count_c = arg.count("c")

    # Tm (°C ) = 2 х (A+T) + 4 х (G+C)
    result = (2 * (count_a + count_t)) + (4 * (count_g + count_c))
    return result


# is_primer — является ли поданная последовательность праймером
def is_primer(arg) -> bool:
    """
    Проверяет, соответствует ли поданная последовательность параметрам
    праймера: GC-составу > 0.4 и < 0.6, Tm (°C ) > 55°C

    Parameters
    ----------
    arg : str
        нуклеотидная последовательность

    Returns
    -------
    bool
        True если является праймером, False если не является
    """
    arg = arg.lower()

    if len(arg) < 16 or len(arg) > 30 or arg[-1].lower() not in ["g", "c"]:
        return False
    else:
        count_g = arg.count("g")
        count_c = arg.count("c")

        temp = annealing_temperature(arg)
        gc = (count_g + count_c) / len(arg)

    return 0.4 < gc < 0.6 or 55 < temp


# is_nucleotide — является ли поданная последовательность нуклеотидной
def is_nucleotide(*args) -> bool:
    """
    Проверяет, является ли поданное значение ДНК или РНК

    Parameters
    ----------
    *args
        список из последовательностей

    Returns
    -------
    bool
        True если является ДНК или РНК, False если не является
    """
    for arg in args:
        if "u" in arg.lower() and "t" in arg.lower():
            return False
        for letter in arg:
            if letter.lower() not in "atgcu":
                return False
    return True
