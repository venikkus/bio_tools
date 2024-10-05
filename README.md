- [russian version](#russian-version)
- [english version](#english-version)

# russian-version
---
# Bio tools

**Bio tools** это набор пакетов для фильтрации FASTQ файлов и работы
с ДНК или РНК последовательностями.

---
Авторы:
* **Разработка:** [*venikkus*](https://github.com/venikkus),
  Санкт-Петербург, Россия.
* **Идея, руководители:**,
  [Институт Биоинформатики](https://bioinf.me/), Санкт-Петербург,
  Россия.
---
## Content

* [Установка](#Установка)
* [Инструкция по использованию](#Инструкция-по-использованию)
* [Контакты](#Контакты)
* [Функции](#Функции)

---
## Установка

Вам нужно скачать папку на ваш ПК. Вы можете скачать её напрямую или
клонировать репозиторий на ваш компьютер.

```bash
git clone git@github.com:venikkus/bio_tools.git
```

## Инструкция по использованию
Для работы с пакетом вы должны импортировать основной скрипт
(bio_tools.py) и вызвать любую из функций


## Функции

Программа содержит две главные функции и дополнительные модули к ним.
1. **run_dna_rna_tools** принимает строку и операцию и выполняет одно из
   действий.
   - transcribe — транскрибирует последовательность
   - reverse — разворачивает последовательность
   - complement — делает из прямой комплементарную последовательность
   - reverse_complement — делает из прямой обратную комплементарную
     последовательность
   - annealing_temperature — выполняет подсчёт температуры отжига
     праймера
   - is_palindrome — проверяет, является ли поданная последовательность
     биопалиндромом
   - is_primer — проверяет, является ли поданная последовательность
     праймером
   - is_nucleotide — проверяет, является ли поданная последовательность
     нуклеотидной


2. **filter_fastq** принимает словарь с ридами и показателями качества
   прочтения и считает для них
   - GC-состав рида.
   - качество рида.


## Контакты
Сообщайте о проблемах напрямую в
[систему отслеживания ошибок](https://github.com/venikkus/bio_tools/issues)
GitHub.<br/>  
Или по почте [niksamusik@gmail.com](mailto:niksamusik@gmail.com)


### english-version
---
# Bio tools

**Bio tools** is a set of packages for filtering FASTQ files and
processing DNA or RNA sequences.

---
Authors:
* **Software:** [*venikkus*](https://github.com/venikkus),
Saint-Petersburg, Russia.

* **Idea, supervisor:**,
  [Bioinformatics Institute](https://bioinf.me/en), Saint-Petersburg,
  Russia.

---
## Content

* [Installation](#installation)
* [Running instructions](#running-instructions)
* [Features](#features)
* [Contact](#contact)

## Installation

You need to download the folder to your PC. You can download directly or
clone the repository to your computer.

```bash
git clone git@github.com:venikkus/bio_tools.git
```
---
## Running instructions

To work with the package, you can import the main script (bio_tools)
and call any of the functions.

## Features

The program contains two main functions and additional modules for them.
1. **run_dna_rna_tools** takes a string and an operation and performs
   one of actions.

    - transcribe — transcribes a sequence
    - complement — converts a forward sequence to a complementary
    - reverse — reverses a sequence
    -  reverse_complement — converts a forward sequence to a reverse
       complementary sequence
    -  annealing_temperature — calculates the annealing temperature of a
       sequence
    -  is_palindrome — checks the supplied sequence is a biopalindrome
    -  is_nucleotide — checks the supplied sequence is a nucleotide
    -  is_primer — checks the supplied sequence is a primer

2. **filter_fastq** takes a dictionary with reads and reading quality
   indicators and calculates for them:

    - the GC composition of a read.
    - the quality of a read.


## Contact

Please report any problems directly to the GitHub
[issue tracker](https://github.com/venikkus/bio_tools/issues).<br/>
Also, you can send your feedback to
[niksamusik@gmail.com](mailto:niksamusik@gmail.com).

### english-version
---
# Bio tools

**Bio tools** is a set of packages for filtering FASTQ files and
processing DNA or RNA sequences.

---
Authors:
* **Software:** [*venikkus*](https://github.com/venikkus), Saint-Petersburg, Russia.
* **Idea, supervisor:** *, [Bioinformatics Institute](https://bioinf.me/en),
  Saint-Petersburg, Russia.

---
## Content

* [Installation](#installation)
* [Running instructions](#running-instructions)
* [Features](#features)
* [Contact](#contact)

## Installation

You need to download the folder to your PC. You can download directly or
clone the repository to your computer.

```bash
git clone git@github.com:venikkus/bio_tools.git
```
---
## Running instructions

To work with the package, you can import the main script (bio_tools)
and call any of the functions.

## Features

The program contains two main functions and additional modules for them.
1. run_dna_rna_tools takes a string and an operation and performs one of
   actions.

    - transcribe — transcribes a sequence
    - complement — converts a forward sequence to a complementary
    - reverse — reverses a sequence
    -  reverse_complement — converts a forward sequence to a reverse
       complementary sequence
    -  annealing_temperature — calculates the annealing temperature of a
       sequence
    -  is_palindrome — checks the supplied sequence is a biopalindrome
    -  is_nucleotide — checks the supplied sequence is a nucleotide
    -  is_primer — checks the supplied sequence is a primer

2. filter_fastq takes a dictionary with reads and reading quality
   indicators and calculates for them:

    - the GC composition of a read.
    - the quality of a read.


## Contact

Please report any problems directly to the GitHub
[issue tracker](https://github.com/venikkus/bio_tools/issues).<br/>
Also, you can send your feedback to
[niksamusik@gmail.com](mailto:niksamusik@gmail.com).
