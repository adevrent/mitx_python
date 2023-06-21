# Hangman game
#

# -----------------------------------
# Helper code
# You don't need to understand this helper code,
# but you will have to know how to use the functions
# (so be sure to read the docstrings!)

import random
import string

WORDLIST_FILENAME = "Unit 3/words.txt"

# remaining_list = [letter for letter in 'abcdefghijklmnopqrstuvwxyz']
def loadWords():
    """
    Returns a list of valid words. Words are strings of lowercase letters.
    
    Depending on the size of the word list, this function may
    take a while to finish.
    """
    print("Loading word list from file...")
    # inFile: file
    inFile = open(WORDLIST_FILENAME, 'r')
    # line: string
    line = inFile.readline()
    # wordlist: list of strings
    wordlist = line.split()
    print("  ", len(wordlist), "words loaded.")
    return wordlist

def chooseWord(wordlist):
    """
    wordlist (list): list of words (strings)

    Returns a word from wordlist at random
    """
    return random.choice(wordlist)

# end of helper code
# -----------------------------------

# Load the list of words into the variable wordlist
# so that it can be accessed from anywhere in the program
wordlist = loadWords()

def isWordGuessed(secretWord, lettersGuessed):
    '''
    secretWord: string, the word the user is guessing
    lettersGuessed: list, what letters have been guessed so far
    returns: boolean, True if all the letters of secretWord are in lettersGuessed;
      False otherwise
    '''
    compiled_word = []
    for letter in secretWord:
      if letter in lettersGuessed:
        compiled_word.append(letter)
    compiled_word = "".join(compiled_word)
    return compiled_word == secretWord


def getGuessedWord(secretWord, lettersGuessed):
    '''
    secretWord: string, the word the user is guessing
    lettersGuessed: list, what letters have been guessed so far
    returns: string, comprised of letters and underscores that represents
      what letters in secretWord have been guessed so far.
    '''
    partial_word = []
    for letter in secretWord:
      if letter in lettersGuessed:
        partial_word.append(letter)
      else:
        partial_word.append("_")
    return "".join(partial_word)



def getAvailableLetters(lettersGuessed):
    '''
    lettersGuessed: list, what letters have been guessed so far
    returns: string, comprised of letters that represents what letters have not
      yet been guessed.
    '''
    remaining_list = [letter for letter in string.ascii_lowercase]
    for letter in lettersGuessed:
      if letter in remaining_list:
        remaining_list.remove(letter)
    return "".join(remaining_list)

def hangman(secretWord):
    '''
    secretWord: string, the secret word to guess.

    Starts up an interactive game of Hangman.

    * At the start of the game, let the user know how many 
      letters the secretWord contains.

    * Ask the user to supply one guess (i.e. letter) per round.

    * The user should receive feedback immediately after each guess 
      about whether their guess appears in the computers word.

    * After each round, you should also display to the user the 
      partially guessed word so far, as well as letters that the 
      user has not yet guessed.

    Follows the other limitations detailed in the problem write-up.
    '''
    guess_left = 8
    lettersGuessed = []
    print("Welcome to the game, Hangman!")
    print("I am thinking of a word that is", len(secretWord), "letters long.")
    print("-----------")
    
    while guess_left > 0:
      print("You have", guess_left, "guesses left")
      print("Available Letters:", getAvailableLetters(lettersGuessed))
      guess_input = (input("Please guess a letter: ")).lower()
      if guess_input in lettersGuessed:
        print("Oops! You've already guessed that letter: ", getGuessedWord(secretWord, lettersGuessed))
        print("-----------")
      elif guess_input in secretWord:
        lettersGuessed.append(guess_input)
        print("Good guess: ", getGuessedWord(secretWord, lettersGuessed))
        print("-----------")
      else:
        guess_left = guess_left - 1
        lettersGuessed.append(guess_input)
        print("Oops! That letter is not in my word: ", getGuessedWord(secretWord, lettersGuessed))
        print("-----------")
      if isWordGuessed(secretWord, lettersGuessed):
        break
    
    if isWordGuessed(secretWord, lettersGuessed):
      print("Congratulations, you won!")
    else:
      print("Sorry, you ran out of guesses. The word was", secretWord)
    



# When you've completed your hangman function, uncomment these two lines
# and run this file to test! (hint: you might want to pick your own
# secretWord while you're testing)

secretWord = chooseWord(wordlist).lower()
#secretWord = "c"
hangman(secretWord)
