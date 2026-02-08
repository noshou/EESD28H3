library(shiny)

#' Asymmetric Lotka-Volterra Competition Model
#'
#' Numerically solves the Lotka-Volterra competition equations using Euler's method
#' for two competing species with asymmetric competition coefficients.
#'
#' Equations:
#'    dN1/dt = r1 * N1 * (1 - (N1 + alpha_1*N2) / k1)
#'    dN2/dt = r2 * N2 * (1 - (N2 + alpha_2*N1) / k2)
#'
#' @param dt Numeric. Time step for Euler integration (e.g., 0.05)
#' @param tm Numeric Maximum time (goes from 0->tm w/ steps by dt)
#' @param r1 Numeric. Intrinsic growth rate of species 1
#' @param r2 Numeric. Intrinsic growth rate of species 2
#' @param k1 Numeric. Carrying capacity of species 1
#' @param k2 Numeric. Carrying capacity of species 2
#' @param n0_1 Numeric. Initial population size of species 1
#' @param n0_2 Numeric. Initial population size of species 2
#' @param alpha12 Numeric. Competition coefficient (effect of sp1 on itself relative to sp2)
#' @param alpha21 Numeric. Competition coefficient (effect of sp2 on itself relative to sp1)
#' @param amp_alpha_12 Numeric. Amplitude variation of alpha 1
#' @param amp_alpha_21 Numeric. Amplitude variation of alpha 2
#' @param phs_alpha_12 Numeric. Phase shift variation of alpha 1 (can be negative or positive)
#' @param phs_alpha_21 Numeric. Phase shift variation of alpha 2
#' @param prd_alpha_12 Numeric. Period variation of alpha 1
#' @param prd_alpha_21 Numeric. Period variation of alpha 2
#' @return Data frame with three columns:
#'    Time - time points from input vector t
#'    N1 - population size of species 1 at each time point
#'    N2 - population size of species 2 at each time point
N_asym <- function(
	dt, 
	tm,
	r1, 
	r2, 
	k1, 
	k2, 
	n0_1, 
	n0_2, 
	alpha12, 
	alpha21,
	amp_alpha_12,
	amp_alpha_21, 
	phs_alpha_12, 
	phs_alpha_21,
	prd_alpha_12,
	prd_alpha_21
) {

	# initialize time span
	t <- seq(from=0, to=tm, by=dt)

	# initialize vectors
	rows <- length(t)
	n1_vec <- rep(0, rows)
	n1_vec[1] <- n0_1
	n2_vec <- rep(0, rows)
	n2_vec[1] <- n0_2

	# floating point error control; 
	# if prd = 2*pi prd == 1 (ie no period shift)
	if (abs(prd_alpha_12 - 2*pi) <= 1.0e-06) {
		omega_1 <- 1
	} 
	else if (prd_alpha_12 == 0) {
		omega_1 <- 0
	}
	else { # convert to angular freq.
		omega_1 <- 2*pi/prd_alpha_12
	}
	if (abs(prd_alpha_21 - 2*pi) <= 1.0e-06) {
		omega_2 <- 1
	} 
	else if (prd_alpha_21 == 0) {
		omega_2 <- 0
	}
	else { # convert to angular freq
		omega_2 <- 2*pi/prd_alpha_21
	}

	# create time varying alpha vectors
	alpha12_var <- alpha12 * (1 + amp_alpha_12 * sin(omega_1*(t+phs_alpha_12)))
	alpha21_var <- alpha21 * (1 + amp_alpha_21 * sin(omega_2*(t+phs_alpha_21)))

	# calculate n1 and n2
	for (i in 2:rows) {

		# load previous values 
		n1_prev <- n1_vec[i-1]
		n2_prev <- n2_vec[i-1]
		a1_prev <- alpha12_var[i-1]
		a2_prev <- alpha21_var[i-1] 
		
		n1_vec[i] <- n1_prev + dt * r1 * n1_prev * ((k1 - (a1_prev * n1_prev + n2_prev))/k1)
		n2_vec[i] <- n2_prev + dt * r2 * n2_prev * ((k2 - (n1_prev + a2_prev * n2_prev))/k2)
	}
	
	# output as df
	df <- data.frame(Time = t, N1 = n1_vec, N2 = n2_vec)
	return(df)
}

# Frontend
ui <- fluidPage(
	withMathJax(), # enable latex formatting
	
	titlePanel("Lotka-Volterra Competition Model with Variation"),
	
	sidebarLayout(
		sidebarPanel(
			h4("Preset Scenarios"),
			actionButton("scenario_a", "Scenario A: Species 1 wins"),
			actionButton("scenario_b", "Scenario B: Species 2 wins"),
			actionButton("scenario_c", "Scenario C: Coexistence"),
			actionButton("scenario_d", "Scenario D: Unstable"),
			hr(),

			h4("Initial Population Sizes"),
			numericInput("n1_0", "\\(N_{1,0}\\):", value = 100, min = 0),
			numericInput("n2_0", "\\(N_{2,0}\\):", value = 100, min = 0),
			hr(),

			h4("Carrying Capacities"),
			numericInput("k1", "\\(K_1\\):", value = 500, min = 0),
			numericInput("k2", "\\(K_2\\):", value = 500, min = 0),
			hr(),

			h4("Competition Coefficients (\\(\\alpha\\))"),
			numericInput("alpha12", "\\(\\alpha_{1,2}\\):", value = 1.1, min = 0, step = 0.1),
			numericInput("alpha21", "\\(\\alpha_{2,1}\\):", value = 0.9, min = 0, step = 0.1),
			hr(),

			h4("Growth rates (\\(r\\))"),
			numericInput("r1", "\\(r_1\\):", value = 0.5, min = 0, step = 0.1),
			numericInput("r2", "\\(r_2\\):", value = 0.5, min = 0, step = 0.1),
			hr(),

			h4("Variation"),
			numericInput("amp_alpha_12", "Amplitude of \\(\\alpha_{1,2}\\):", value = 1, min = 0, step = 0.05),
			numericInput("amp_alpha_21", "Amplitude of \\(\\alpha_{2,1}\\):", value = 1, min = 0, step = 0.05),
			numericInput("prd_alpha_12", "Period of \\(\\alpha_{2,1}\\):", value = 2*pi, min = 0, step = 1),
			numericInput("prd_alpha_21", "Period of \\(\\alpha_{2,1}\\):", value = 2*pi, min = 0, step = 1),
			numericInput("phs_alpha_12", "Phase Shift of \\(\\alpha_{1,1}\\):", value = 0, min=-50, max=50, step = 0.5),
			numericInput("phs_alpha_21", "Phase Shift of \\(\\alpha_{2,1}\\):", value = 0, min=-50, max=50, step = 0.5),
			hr(),

			h4("Simulation Settings"),
			actionButton("exec", "Run Simulation", class = "btn-primary"),
			actionButton("scenario_default","Reset Parameters")
		),
		
		# Populations, and Isoclines
		mainPanel(
			tabsetPanel(
				
				tabPanel("Populations",
						wellPanel(
						h4("Display Options"),
							fluidRow(
								column(4, checkboxInput("show_N1", "Show \\(N_1\\)", value = TRUE)),
								column(4, checkboxInput("show_N2", "Show \\(N_2\\)", value = TRUE)),
								column(4, checkboxInput("show_K", "Show Carrying Capacities", value = TRUE))
							)
						),
						plotOutput("timePlot", height = "500px"),
						hr(),
						h4("Final Population Sizes:"),
						uiOutput("finalPops")),
				
				tabPanel("Isoclines",
						wellPanel(
						h4("Display Options"),
							fluidRow(
								column(4, checkboxInput("show_K_iso", "Show Carrying Capacities", value = TRUE))
							)
						),
						plotOutput("isoclinePlot", height = "500px"),
						hr(),
						helpText("Nullclines (isoclines) show where each population has zero growth. 
								Intersection points represent equilibria."))
			)
		)
	)
)

# Backend

server <- function(input, output, session) {
	# Default Scenario
	observeEvent(input$scenario_default, {
		updateNumericInput(session, "k1", value=500)
		updateNumericInput(session, "k2", value=500)
		updateNumericInput(session, "alpha12", value=1.1)
		updateNumericInput(session, "alpha21", value=0.9)
		updateNumericInput(session, "r1", value=0.5)
		updateNumericInput(session, "r2", value=0.5)
		updateNumericInput(session, "n1_0", value=100)
		updateNumericInput(session, "n2_0", value=100)
		updateNumericInput(session, "amp_alpha_12", value=1)
		updateNumericInput(session, "amp_alpha_21", value=1)
	})

	# Scenario A
	observeEvent(input$scenario_a, {
		updateNumericInput(session, "k1", value=400)
		updateNumericInput(session, "k2", value=250)
		updateNumericInput(session, "alpha12", value=0.5)
		updateNumericInput(session, "alpha21", value=2)
		updateNumericInput(session, "r1", value=0.5)
		updateNumericInput(session, "r2", value=0.5)
		updateNumericInput(session, "n1_0", value=50)
		updateNumericInput(session, "n2_0", value=50)
		updateNumericInput(session, "amp_alpha_12", value=0)
		updateNumericInput(session, "amp_alpha_21", value=0)
	})

	# Scenario B
	observeEvent(input$scenario_b, {
		updateNumericInput(session, "k1", value=400)
		updateNumericInput(session, "k2", value=250)
		updateNumericInput(session, "alpha12", value=2)
		updateNumericInput(session, "alpha21", value=0.5)
		updateNumericInput(session, "r1", value=0.5)
		updateNumericInput(session, "r2", value=0.5)
		updateNumericInput(session, "n1_0", value=50)
		updateNumericInput(session, "n2_0", value=50)
		updateNumericInput(session, "amp_alpha_12", value=0)
		updateNumericInput(session, "amp_alpha_21", value=0)
	})

	# Scenario C
	observeEvent(input$scenario_c, {
		updateNumericInput(session, "k1", value=400)
		updateNumericInput(session, "k2", value=250)
		updateNumericInput(session, "alpha12", value=2.5)
		updateNumericInput(session, "alpha21", value=1.5)
		updateNumericInput(session, "r1", value=0.5)
		updateNumericInput(session, "r2", value=0.5)
		updateNumericInput(session, "n1_0", value=50)
		updateNumericInput(session, "n2_0", value=50)
		updateNumericInput(session, "amp_alpha_12", value=0)
		updateNumericInput(session, "amp_alpha_21", value=0)
	})

	# Scenario D
	observeEvent(input$scenario_d, {
		updateNumericInput(session, "k1", value=400)
		updateNumericInput(session, "k2", value=250)
		updateNumericInput(session, "alpha12", value=0.5)
		updateNumericInput(session, "alpha21", value=0.5)
		updateNumericInput(session, "r1", value=0.5)
		updateNumericInput(session, "r2", value=0.5)
		updateNumericInput(session, "n1_0", value=50)
		updateNumericInput(session, "n2_0", value=50)
		updateNumericInput(session, "amp_alpha_12", value=0)
		updateNumericInput(session, "amp_alpha_21", value=0)
	})
	# Simulation
	sim <- eventReactive(input$exec, {
		df <- N_asym(
			dt 			= 0.05,
			tm 			= 100,
			r1 			= input$r1, 
			r2 			= input$r2, 
			k1 			= input$k1,
			k2 			= input$k2,
			n0_1 		= input$n1_0,
			n0_2 		= input$n2_0,
			alpha12     = input$alpha12,
			alpha21     = input$alpha21,
			amp_alpha_12 = input$amp_alpha_12,
			amp_alpha_21 = input$amp_alpha_21,
			phs_alpha_12 = input$phs_alpha_12, 
			phs_alpha_21 = input$phs_alpha_21,
			prd_alpha_12 = input$prd_alpha_12, 
			prd_alpha_21 = input$prd_alpha_21
		)
		return(df)
	})
	
	output$timePlot <- renderPlot({
		req(sim())		
		data <- sim()
		
		# Determine y-axis limits
		y_max <- 0
		if(input$show_N1) y_max <- max(y_max, max(data$N1))
		if(input$show_N2) y_max <- max(y_max, max(data$N2))
		if(input$show_K) y_max <- max(y_max, input$k1, input$k2)
		if(y_max == 0) y_max <- 100  # Default if nothing selected
		
		# Create plot
		plot(data$Time, data$N1, 
			type = ifelse(input$show_N1, "l", "n"),  # "n" = no plotting if unchecked
			col = "blue", lwd = 2,
			 ylim = c(0, y_max * 1.1),
			xlab = "Time", 
			ylab = "Population Size",
			main = "Lotka-Volterra Competition Dynamics")
		
		# Add N2 if selected
		if(input$show_N2) {
			lines(data$Time, data$N2, col = "red", lwd = 2)
		}
		
		# Add carrying capacities if selected
		if(input$show_K) {
			abline(h = input$k1, col = "lightblue", lwd = 2, lty = 3)
			abline(h = input$k2, col = "pink", lwd = 2, lty = 3)
		}
		
		# Build legend
		legend_items <- c()
		legend_cols <- c()
		legend_lty <- c()
		
		if(input$show_N1) {
			legend_items <- c(legend_items, "Species 1")
			legend_cols <- c(legend_cols, "blue")
			legend_lty <- c(legend_lty, 1)
		}
		if(input$show_N2) {
			legend_items <- c(legend_items, "Species 2")
			legend_cols <- c(legend_cols, "red")
			legend_lty <- c(legend_lty, 1)
		}
		if(input$show_K) {
			legend_items <- c(legend_items, expression(K[1]), expression(K[2]))
			legend_cols <- c(legend_cols, "lightblue", "pink")
			legend_lty <- c(legend_lty, 3, 3)
		}
		
		if(length(legend_items) > 0) {
			legend("right", 
					legend = legend_items, 
					col = legend_cols, 
					lty = legend_lty,
					lwd = 2, 
					cex = 0.8,
					bty = "n")
		}
	})
	
output$isoclinePlot <- renderPlot({
	req(sim())		
	data <- sim()
	
	# Calculate plot limits
	xlim <- c(0, max(input$k1, input$k2/input$alpha21) * 1.15)
	ylim <- c(0, max(input$k2, input$k1/input$alpha12) * 1.15)
	
	# Create plot with proper axes
	plot(c(input$k1, 0), c(0, input$k1/input$alpha12),
		type = "l", col = "blue", lwd = 2,
		xlim = xlim, ylim = ylim,
		xlab = expression(N[1]),
		ylab = expression(N[2]),
		main = "Isoclines (Zero Growth Lines)",
		xaxt = "n", yaxt = "n")
	
	grid(col = "lightgray", lty = "dotted")
	
	# Add Species 2 isocline
	lines(c(0, input$k2/input$alpha21), c(input$k2, 0), col = "red", lwd = 2)
	
	# Add trajectory
	lines(data$N1, data$N2, col = "darkgreen", lwd = 2)
	points(data$N1[1], data$N2[1], pch = 21, bg = "green", cex = 2)
	points(data$N1[nrow(data)], data$N2[nrow(data)], pch = 22, bg = "red", cex = 2)
	
	# Add custom axes if checkbox selected
	if (input$show_K_iso) {
		axis(1, at = c(input$k2/input$alpha21, input$k1), 
			labels = c(expression(K[2]/alpha[21]), expression(K[1])),
			las = 1)
		axis(2, at = c(input$k1/input$alpha12, input$k2), 
			labels = c(expression(K[1]/alpha[12]), expression(K[2])),
			las = 2)
	} else {
		axis(1)
		axis(2)
	}
	
	# Legend
	legend("topright",
			legend = c(expression(paste(N[1], " isocline")),
					expression(paste(N[2], " isocline")),
					"Trajectory", "Start", "End"),
			col = c("blue", "red", "darkgreen", "darkgreen", "darkred"),
			lty = c(1, 1, 1, NA, NA),
			pch = c(NA, NA, NA, 21, 22),
			pt.bg = c(NA, NA, NA, "green", "red"),
			lwd = 2,
			bty = "n")
})	
	output$finalPops <- renderUI({
		req(sim())		
		data <- sim()
		withMathJax(
			HTML(paste0(
				"\\(N_1(t_{final})\\) = ", round(data$N1[nrow(data)], 2), "<br/>",
				"\\(N_2(t_{final})\\) = ", round(data$N2[nrow(data)], 2)
			))
		)
	})
}

shinyApp(ui, server, options = list(launch.browser = TRUE))