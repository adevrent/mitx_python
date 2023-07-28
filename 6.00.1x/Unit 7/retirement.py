class retirementFund(object):
    def __init__(self, m, r, r_a, c_a) -> None:
        self.m = m  # monthly amount saved in $.
        self.r = r  # percentage earned from investments each month (decimal.)
        self.r_a = r_a  # retirement age (in years.)
        self.c_a = c_a  # current age (in years.)
        self.total_savings = 0  # Total savings at retirement (in $.)
        self.monthly_savings = []  # Total accumulated savings in each month.
        
    def calculate_savings(self):
        """Calculates total savings at
        the age of retirement.
        Also calculates total accumulated 
        savings in each month.

        Args:
            m (float): monthly amount saved in $.
            r (float): percentage earned from investments each month (decimal.)
            r_a (int): retirement age (in years.)
            c_a (int): current age (in years.)
        """
        total_savings = 0
        monthly_savings = []
        for month in range((self.r_a - self.c_a) * 12):
            total_savings *= (1 + self.r)
            total_savings += self.m
            monthly_savings.append(total_savings)
        self.total_savings = total_savings
        self.monthly_savings = monthly_savings
        
        return total_savings
    
    def getMonthlySavings(self):
        return self.monthly_savings
    
    def setRetirementAge(self, r_a):
        self.r_a = r_a
        
    def setCurrentAge(self, c_a):
        self.c_a = c_a
    
    def setMonthlySaving(self, m):
        self.m = m
    
    def setMonthlyRatio(self, r):
        self.r = r
    
    def __str__(self) -> str:
        return "By saving " + "$" + str(self.m) + " per month, at age " + str(self.r_a) + " you are going to have " + "$" + str(round(self.total_savings, 2))
    
    def getInfo(self):
        print("Current monthly savings:", "$" + str(self.m))
        print("Current monthly return rate:",  str(100 * self.r) + "%")
        print("Current age:", str(self.c_a))
        print("Current retirement age:", str(self.r_a))
        print("Current expected savings at retirement:", "$" + str(self.total_savings))
        