# Candy Distribution Algorithm
def distribute_candies(ratings):
    n = len(ratings)
    candies = [1] * n  # Initialize with 1 candy each
    
    # Forward pass: Compare with left neighbor
    for i in range(1, n):
        if ratings[i] > ratings[i-1]:
            candies[i] = candies[i-1] + 1
    
    # Backward pass: Compare with right neighbor
    for i in range(n-2, -1, -1):
        if ratings[i] > ratings[i+1]:
            candies[i] = max(candies[i], candies[i+1] + 1)
    
    return candies, sum(candies)

# Example with 6 random numbers
ratings = [5, 2, 8, 4, 1, 6]
candies, total = distribute_candies(ratings)

print("Ratings:", ratings)
print("Candies:", candies)
print("Total candies needed:", total)