from modules import dna_rna_tools_modules, filter_fastq_modules

"""
Модуль предоставляет функции для обработки
ДНК/РНК последовательностей и фильтрации данных FASTQ.

Функции:
- run_dna_rna_tools: Выполняет операции над последовательностями ДНК/РНК.
- filter_fastq: Фильтрует FASTQ файлы по заданным
параметрам качества, длины и GC-содержания.

Автор: Ника Самусик
"""


def run_dna_rna_tools(*args):
    """
    Выполняет операцию с последовательностями ДНК или РНК
    из модуля dna_rna_tools_modules

    Parameters
    ----------
    *args : tuple, str
        последовательности ДНК/РНК и элемент с названием операции.
        Операция должна быть передана последним аргументом.

    Returns
    -------
    answer : list, str
        Результаты или результат выполнения функций из модуля
        dna_rna_tools_modules. Предупреждение о неправильном вводе,
        если поданные данные не содержат операции или подают неверную операцию

    Raises
    -------
    ValueError
        если последовательности не являются ДНК/РНК или
    KeyError
        если операция не поддерживается
    """
    if len(args) == 1:
        raise ValueError("There is no operation or sequence")

    operation = args[-1]
    args = args[:-1]

    # check if it is a DNA or RNA at all
    if not dna_rna_tools_modules.is_nucleotide(*args):
        raise ValueError("This is not a DNA/RNA sequence at all")

    if operation in [
        "transcribe",
        "reverse",
        "complement",
        "reverse_complement",
        "annealing_temperature",
        "is_palindrome",
        "is_primer",
    ]:
        answer = [
            getattr(dna_rna_tools_modules, operation)(arg) for arg in args
        ]
        return answer[0] if len(answer) == 1 else answer
    else:
        raise KeyError(f"Operation {operation} is not supported.")


def filter_fastq(
    fastq_file,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0
):
    """
    Фильтрует последовательности на соответствие параметрам качества,
    длины и GC-состава.

    Parameters
    ----------
    fastq_file : dict
        входные данные: словарь, где ключ — имя последовательности,
        значение — кортеж (строка последовательности, строка качества).
    gc_bounds : int, tuple, default: 0, 100
        пороги содержания GC в последовательности.
        Можно задать только верхний порог
    length_bounds : int, tuple, default: 0, 2**32
        пороги длины последовательности. Можно задать только верхний порог
    quality_threshold : int, default: 0
        Минимальный порог качества для фильтрации последовательностей

    Returns
    -------
    dict
        Отфильтрованные данные FASTQ

    Raises
    -------
    ValueError
        если значения GC-состава или длины выходят за пределы.
    """
    filtered_fastq = {}
    for seq_name, seq_data in fastq_file.items():
        seq, quality = seq_data

        gc_content = filter_fastq_modules.calculate_gc_bounds(seq)
        seq_length = len(seq)
        q_score = filter_fastq_modules.calculate_quality_threshold(quality)

        if isinstance(gc_bounds, int):
            lower_gc_threshold, upper_gc_threshold = 0, gc_bounds
        else:
            lower_gc_threshold, upper_gc_threshold = gc_bounds

        if lower_gc_threshold < 0 or upper_gc_threshold > 100:
            raise ValueError("GC composition limits must be between 0 and 100")

        if isinstance(length_bounds, int):
            lower_length_threshold, upper_length_threshold = 0, length_bounds
        else:
            lower_length_threshold, upper_length_threshold = length_bounds

        if (
            lower_length_threshold < 0
            or upper_length_threshold < lower_length_threshold
        ):
            raise ValueError("Length limits are specified incorrectly.")

        if (
            lower_gc_threshold <= gc_content <= upper_gc_threshold
            and lower_length_threshold <= seq_length <= upper_length_threshold
            and q_score >= quality_threshold
        ):
            filtered_fastq[seq_name] = (seq, quality)

    return filtered_fastq
