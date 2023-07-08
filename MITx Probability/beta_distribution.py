import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, uniform
import ipywidgets as widgets
from IPython.display import display, clear_output

def calculate_posterior_pdf(theta_prior, prior_type, N, K):
    # Generate theta values between 0 and 1
    theta_values = np.linspace(0, 1, 1000)

    # Calculate the likelihood function
    likelihood = theta_values**K * (1 - theta_values)**(N - K)

    # Calculate the posterior using Bayes' rule
    if prior_type == 'beta':
        posterior = likelihood * theta_prior.pdf(theta_values)
    elif prior_type == 'uniform':
        posterior = likelihood * 1  # Uniform prior PDF is constant (1)
    else:
        raise ValueError('Invalid prior type. Choose either "beta" or "uniform".')

    # Normalize the posterior PDF
    posterior /= np.trapz(posterior, theta_values)  # Normalize by dividing by the integral

    return theta_values, posterior

def update_plot(prior_type, N, alpha, beta_value):
    # Simulated data
    K = np.random.binomial(N, true_bias)  # Number of observed heads

    # Calculate the posterior PDF
    if prior_type == 'beta':
        theta_prior = beta(alpha, beta_value)
    elif prior_type == 'uniform':
        theta_prior = uniform()
    else:
        raise ValueError('Invalid prior type. Choose either "beta" or "uniform".')

    theta_values, posterior = calculate_posterior_pdf(theta_prior, prior_type, N, K)

    # Clear the previous plot
    clear_output(wait=True)

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(theta_values, theta_prior.pdf(theta_values), 'r-', label='Prior')
    plt.plot(theta_values, posterior, 'b-', label='Posterior')
    plt.xlabel('Theta')
    plt.ylabel('Probability Density')
    plt.title('Prior and Posterior Distributions of Theta')
    plt.legend()
    plt.show()

    # Re-display the widgets
    display(widgets_box)

# True bias of the coin
true_bias = 0.3

# Create the interactive widgets
prior_type_widget = widgets.Dropdown(options=['uniform', 'beta'], value='uniform', description='Prior Type:') 
N_widget = widgets.IntSlider(min=1, max=100, value=50, description='N:', continuous_update=False)
alpha_widget = widgets.FloatSlider(min=0.01, max=10, value=1, description='Alpha:', continuous_update=False)
beta_value_widget = widgets.FloatSlider(min=0.01, max=10, value=1, description='Beta:', continuous_update=False)

# Define the function to be called when the widgets change
def on_widgets_change(change):
    prior_type = prior_type_widget.value
    N = N_widget.value
    alpha = alpha_widget.value
    beta_value = beta_value_widget.value
    update_plot(prior_type, N, alpha, beta_value)

# Register the widget event handlers
prior_type_widget.observe(on_widgets_change, 'value') 
N_widget.observe(on_widgets_change, 'value')
alpha_widget.observe(on_widgets_change, 'value')
beta_value_widget.observe(on_widgets_change, 'value')

# Create a box to hold the widgets
widgets_box = widgets.VBox([prior_type_widget, N_widget, alpha_widget, beta_value_widget])

# Display the widgets
display(widgets_box)

# Initialize the plot
update_plot(prior_type_widget.value, N_widget.value, alpha_widget.value, beta_value_widget.value)