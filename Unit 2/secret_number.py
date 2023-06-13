print("Please think of a number between 0 and 100!")
low = 0
high = 100

while True:
    guess = (low + high) // 2
    print("Is your secret number " + str(guess) + "?")
    answer = input("Enter 'h' to indicate the guess is too high. Enter 'l' to indicate the guess is too low. Enter 'c' to indicate I guessed correctly. ")
    if answer == 'c':
        break
    elif answer == 'h':
        high = guess
    elif answer == 'l':
        low = guess
    else:
        print("Sorry, I did not understand your input.")
print("Game over. Your secret number was:", guess)