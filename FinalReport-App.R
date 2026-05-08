# =============================================================================
# Final Portfolio Shiny App
# Four models: Lotka-Volterra Competition, Eutrophication, Streeter-Phelps,
#              1D Mass Transport (Advection-Diffusion)
# =============================================================================

library(shiny)
library(bslib)
library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

# =============================================================================
# MODEL 1: Lotka-Volterra Competition with time-varying alpha (original code)
# =============================================================================

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

# =============================================================================
# MODEL 2: Eutrophication (logistic growth + Type-III grazing)
# dP/dt = r*P*(1 - P/K) - g*P^2 / (h^2 + P^2)
# =============================================================================
eutroph_sim <- function(r, h, g, K, P0, tmax = 400, dt = 0.1) {
	time <- seq(0, tmax, by = dt)
	P <- numeric(length(time)); P[1] <- P0
	for (i in 2:length(time)) {
		growth  <- r * P[i - 1] * (1 - P[i - 1] / K)
		grazing <- (g * P[i - 1]^2) / (h^2 + P[i - 1]^2)
		P[i] <- P[i - 1] + dt * (growth - grazing)
	}
	data.frame(time = time, P = P)
}

eutroph_dPdt <- function(P, r, h, g, K) {
	r * P * (1 - P / K) - (g * P^2) / (h^2 + P^2)
}

# =============================================================================
# MODEL 3: Streeter-Phelps + parameter fitting via optim()
# =============================================================================
streeter_phelps <- function(time, state, parameters) {
	with(as.list(c(state, parameters)), {
		dL <- -k_d * L
		dD <-  k_d * L - k_a * D
		list(c(dL, dD))
	})
}

sp_simulate <- function(L0, D0, k_d, k_a, times) {
	state <- c(L = L0, D = D0)
	parms <- c(k_d = k_d, k_a = k_a)
	as.data.frame(ode(y = state, times = times, func = streeter_phelps, parms = parms))
}

sp_error <- function(params, data, DO_sat) {
	out <- sp_simulate(L0 = params[3], D0 = params[4], k_d = params[1], k_a = params[2],
					   times = data$time)
	out$DO <- DO_sat - out$D
	BOD_err <- sum((out$L  - data$BOD)^2, na.rm = TRUE)
	DO_err  <- sum((out$DO - data$DO )^2, na.rm = TRUE)
	BOD_err + 5 * DO_err
}

# =============================================================================
# MODEL 4: 1D Advection-Diffusion (Method of Lines)
# Two boundary condition options on the downstream side.
# =============================================================================
transport_model <- function(time, state, parameters) {
	with(as.list(parameters), {
		C <- state
		n <- length(C)
		dC_dt <- numeric(n)
		for (i in 1:n) {
			if (i == 1) {
				# Upstream: Dirichlet, fixed at C_upstream
				C[i] <- C_upstream
				dC_dt[i] <- 0
			} else if (i == n) {
				if (bc_down == "neumann") {
					# Zero gradient: outflow boundary
					C[i] <- C[i - 1]
					dC_dt[i] <- 0
				} else {
					# Dirichlet: fixed at C_downstream
					C[i] <- C_downstream
					dC_dt[i] <- 0
				}
			} else {
				advection <- -u * (C[i] - C[i - 1]) / dx
				diffusion <- D * (C[i - 1] - 2 * C[i] + C[i + 1]) / (dx^2)
				dC_dt[i]  <- advection + diffusion
			}
		}
		list(c(dC_dt))
	})
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
	theme = bs_theme(bootswatch = "journal"),
	withMathJax(),
	titlePanel("Environmental Modelling: Final Portfolio"),

	tabsetPanel(
		id = "tabs",

		# ----- About ---------------------------------------------------------
		tabPanel("About",
			fluidRow(column(10, offset = 1,
				h3("Overview"),
				p("This app implements four environmental models covered through the term:"),
				tags$ul(
					tags$li(strong("Competition:"), " coupled Lotka-Volterra dynamics for two species, with optional time-varying competition coefficients."),
					tags$li(strong("Eutrophication:"), " logistic growth with non-linear grazing producing alternative stable states."),
					tags$li(strong("Streeter-Phelps:"), " parameterized BOD/DO river model fit to imported data via ", code("optim()"), "."),
					tags$li(strong("Mass Transport 1D:"), " advection-diffusion PDE solved by the method of lines.")
				),
				p("Each tab provides editable parameters, state variables (where relevant), and labelled
				   plots with explanations of what to look for.")
			))
		),

		# ----- Competition (original UI preserved) ---------------------------
		tabPanel("Competition",
			sidebarLayout(
				sidebarPanel(
					h4("Preset Scenarios"),
					actionButton("scenario_a", "Scenario A: Species 1 wins"),
					actionButton("scenario_b", "Scenario B: Species 2 wins"),
					actionButton("scenario_c", "Scenario C: Coexistence"),
					actionButton("scenario_d", "Scenario D: Unstable"),
					hr(),

					h4("Initial Population Sizes"),
					numericInput("n1_0", "\\(N_{1,0}\\):", value = 550, min = 0),
					numericInput("n2_0", "\\(N_{2,0}\\):", value = 1, min = 0),
					hr(),

					h4("Carrying Capacities"),
					numericInput("k1", "\\(K_1\\):", value = 500, min = 0),
					numericInput("k2", "\\(K_2\\):", value = 600, min = 0),
					hr(),

					h4("Competition Coefficients (\\(\\alpha\\))"),
					numericInput("alpha12", "\\(\\alpha_{12}\\):", value = 1.4, min = 0, step = 0.1),
					numericInput("alpha21", "\\(\\alpha_{21}\\):", value = 1.2, min = 0, step = 0.1),
					hr(),

					h4("Growth rates (\\(r\\))"),
					numericInput("r1", "\\(r_1\\):", value = 1.68, min = 0, step = 0.01),
					numericInput("r2", "\\(r_2\\):", value = 1.5, min = 0, step = 0.01),
					hr(),

					h4("Variation"),
					numericInput("amp_alpha_12", "Amplitude of \\(\\alpha_{12}\\):", value = 1, min = 0, step = 0.025),
					numericInput("amp_alpha_21", "Amplitude of \\(\\alpha_{21}\\):", value = 1, min = 0, step = 0.025),
					numericInput("prd_alpha_12", "Period of \\(\\alpha_{12}\\):", value = 2.130, min = 0, step = 1),
					numericInput("prd_alpha_21", "Period of \\(\\alpha_{21}\\):", value = 0.5, min = 0, step = 1),
					numericInput("phs_alpha_12", "Phase Shift of \\(\\alpha_{12}\\):", value = 0, step = 0.025),
					numericInput("phs_alpha_21", "Phase Shift of \\(\\alpha_{21}\\):", value = -2.5, step = 0.025),
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
								uiOutput("finalPops"),
								helpText("Population sizes through time. Whichever curve survives at large t
										 identifies the winner; if both stabilize above zero, the species coexist.
										 The dashed lines mark each species' carrying capacity.")),
						
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
										Intersection points represent equilibria. The trajectory starts at the
										green dot and ends at the red square; following it shows which equilibrium
										the chosen initial conditions get pulled towards."))
					)
				)
			)
		),

		# ----- Eutrophication ------------------------------------------------
		tabPanel("Eutrophication",
			sidebarLayout(
				sidebarPanel(
					h4("State variable"),
					numericInput("eu_P0", "Initial \\(P_0\\):", value = 1, min = 0, step = 0.5),
					hr(),
					h4("Parameters"),
					numericInput("eu_r", "Growth rate \\(r\\):", value = 0.5, min = 0, step = 0.05),
					numericInput("eu_h", "Half-saturation \\(h\\):", value = 1.5, min = 0.01, step = 0.1),
					numericInput("eu_g", "Max grazing \\(g\\):", value = 1.5, min = 0, step = 0.1),
					numericInput("eu_K", "Carrying capacity \\(K\\):", value = 15, min = 0.1, step = 1),
					hr(),
					numericInput("eu_tmax", "\\(t_{max}\\):", value = 200, min = 10),
					actionButton("eu_run", "Run Simulation", class = "btn-primary")
				),
				mainPanel(
					h4("\\(P(t)\\): trajectory in time"),
					plotOutput("euTimePlot", height = "350px"),
					helpText("Phytoplankton biomass over time, starting from \\(P_0\\). The curve flattens
							  out as it approaches a stable equilibrium. Try several initial conditions
							  (e.g. \\(P_0 = 0.5,\\, 6,\\, 12\\)) to see whether they all settle on the same
							  equilibrium or get pulled to different stable states."),
					hr(),
					h4("Phase portrait: \\(dP/dt\\) vs \\(P\\)"),
					plotOutput("euPhasePlot", height = "350px"),
					helpText("Zeros of \\(dP/dt\\) (where the curve crosses the dashed line) are
							  equilibria. A zero where the curve goes from positive to negative is
							  stable; the opposite is unstable. With strong grazing the system can
							  display two stable equilibria separated by an unstable one — the classic
							  alternative-stable-states regime.")
				)
			)
		),

		# ----- Streeter-Phelps -----------------------------------------------
		tabPanel("Streeter-Phelps",
			sidebarLayout(
				sidebarPanel(
					h4("Data"),
					fileInput("sp_file", "Upload CSV (columns: time, BOD, DO):",
							  accept = c(".csv")),
					numericInput("sp_DO_sat", "\\(DO_{sat}\\) (mg/L):", value = 9.2, step = 0.1),
					hr(),
					h4("Initial guesses"),
					numericInput("sp_L0", "\\(L_0\\):", value = 10, step = 0.5),
					numericInput("sp_D0", "\\(D_0\\):", value = 2, step = 0.5),
					numericInput("sp_kd", "\\(k_d\\):", value = 1, step = 0.1),
					numericInput("sp_ka", "\\(k_a\\):", value = 1, step = 0.1),
					hr(),
					actionButton("sp_fit", "Fit Model", class = "btn-primary"),
					hr(),
					helpText("A built-in default dataset is loaded automatically. Upload a CSV with columns time, BOD, DO to override it.")
				),
				mainPanel(
					h4("Raw data (no model)"),
					plotOutput("spRawPlot", height = "300px"),
					helpText("Use this plot to eyeball reasonable initial guesses before fitting:
							  \\(L_0\\) is the BOD intercept at \\(t=0\\), \\(D_0 = DO_{sat} - DO(0)\\),
							  and \\(k_d, k_a\\) control how fast the curves decay/recover."),
					hr(),
					h4("Fitted parameters"),
					verbatimTextOutput("spFitText"),
					h4("Optimized model vs. data"),
					plotOutput("spFitPlot", height = "350px"),
					helpText("Solid curves are the optimized model; points are the input data.
							  A close match indicates the optimizer found a good local minimum;
							  if the curves miss the points badly, try different initial guesses.")
				)
			)
		),

		# ----- Mass Transport 1D ---------------------------------------------
		tabPanel("Mass Transport 1D",
			sidebarLayout(
				sidebarPanel(
					h4("Parameters"),
					numericInput("mt_u", "Advection velocity \\(u\\) (m/day):", value = 10, step = 1),
					numericInput("mt_D", "Diffusion coefficient \\(D\\) (m²/day):", value = 10, step = 1),
					hr(),
					h4("Boundary conditions"),
					selectInput("mt_bc", "Downstream boundary:",
								choices = c("Neumann (zero gradient / open outflow)" = "neumann",
											"Dirichlet (fixed concentration)"        = "dirichlet")),
					numericInput("mt_C_up",   "Upstream \\(C\\) (Dirichlet):",   value = 0, step = 0.5),
					numericInput("mt_C_down", "Downstream \\(C\\) (if Dirichlet):", value = 0, step = 0.5),
					hr(),
					h4("Initial impulse"),
					numericInput("mt_spike_pos", "Spike position (box index):", value = 2, min = 1, step = 1),
					numericInput("mt_spike_C",   "Spike concentration:",         value = 3, step = 0.5),
					hr(),
					numericInput("mt_tmax", "\\(t_{max}\\) (days):", value = 50, min = 1),
					actionButton("mt_run", "Run Simulation", class = "btn-primary")
				),
				mainPanel(
					h4("Concentration heatmap"),
					plotOutput("mtHeatPlot", height = "500px"),
					helpText("Color encodes concentration; the \\(x\\)-axis is distance downstream
							  and the \\(y\\)-axis is time. Pure advection (\\(D = 0\\)) shows a
							  diagonal stripe of constant width: the spike translates without spreading.
							  Pure diffusion (\\(u = 0\\)) shows a symmetric plume centred on the spike
							  that widens and fades. With both, the plume is carried downstream and
							  spreads simultaneously. The Neumann boundary lets mass leave the domain
							  freely; the Dirichlet boundary clamps the downstream end.")
				)
			)
		)
	)
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

	# ---------- Competition (original server logic preserved) ---------------
	# Default Scenario
	observeEvent(input$scenario_default, {
		updateNumericInput(session, "k1", value=500)
		updateNumericInput(session, "k2", value=600)
		updateNumericInput(session, "alpha12", value=1.4)
		updateNumericInput(session, "alpha21", value=1.2)
		updateNumericInput(session, "r1", value=1.68)
		updateNumericInput(session, "r2", value=1.5)
		updateNumericInput(session, "n1_0", value=550)
		updateNumericInput(session, "n2_0", value=1)
		updateNumericInput(session, "amp_alpha_12", value=1)
		updateNumericInput(session, "amp_alpha_21", value=1)
		updateNumericInput(session, "phs_alpha_12", 0)
		updateNumericInput(session, "phs_alpha_21", -2.5)
		updateNumericInput(session, "prd_alpha_12", 2.13)
		updateNumericInput(session, "prd_alpha_21", 0.5)
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

	# ---------- Eutrophication ----------------------------------------------
	eu_sim <- eventReactive(input$eu_run, {
		eutroph_sim(r = input$eu_r, h = input$eu_h, g = input$eu_g,
					K = input$eu_K, P0 = input$eu_P0, tmax = input$eu_tmax)
	}, ignoreNULL = FALSE)

	output$euTimePlot <- renderPlot({
		d <- eu_sim()
		ggplot(d, aes(x = time, y = P)) +
			geom_line(linewidth = 1.1, color = "darkgreen") +
			labs(title = "Phytoplankton biomass through time",
				 x = "Time", y = "P (biomass)") +
			theme_minimal(base_size = 13)
	})

	output$euPhasePlot <- renderPlot({
		Pmax <- max(input$eu_K * 1.1, 1)
		Pv <- seq(0, Pmax, length.out = 600)
		dPv <- eutroph_dPdt(Pv, input$eu_r, input$eu_h, input$eu_g, input$eu_K)
		# Detect equilibria via sign changes
		sign_change <- which(diff(sign(dPv)) != 0)
		eq <- Pv[sign_change]
		ggplot(data.frame(P = Pv, dP = dPv), aes(x = P, y = dP)) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
			geom_line(linewidth = 1.1, color = "steelblue") +
			geom_point(data = data.frame(P = eq, dP = rep(0, length(eq))),
					   aes(x = P, y = dP), color = "red", size = 3) +
			labs(title = "Phase portrait: dP/dt vs P",
				 subtitle = paste("Equilibria detected:", length(eq)),
				 x = "P", y = "dP/dt") +
			theme_minimal(base_size = 13)
	})

	# ---------- Streeter-Phelps ---------------------------------------------
	# Hardcoded default dataset (equivalent to Streeter-Phelps_data_opt.csv).
	# Users can override by uploading their own CSV with columns: time, BOD, DO.
	default_sp_data <- data.frame(
		time = seq(0, 10, by = 0.5),
		BOD  = c(11.4572, 8.5425, 5.5334, 2.8612, 2.1335, 2.4497, -0.1247,
				 0.5153, 1.1221, -0.1055, -0.1197, 0.1000, 0.8445, -0.2533,
				 -0.1776, -0.1874, 1.1229, 1.1068, 0.5110, 0.1991, 0.3727),
		DO   = c(6.6472, 6.6968, 7.9063, 7.6041, 8.1051, 8.8027, 8.2755,
				 8.7906, 8.6693, 8.9217, 8.2111, 8.5579, 8.9081, 9.4175,
				 9.1029, 9.1675, 9.3818, 8.9178, 9.2715, 8.9484, 8.6743)
	)

	sp_data <- reactive({
		if (is.null(input$sp_file)) {
			default_sp_data
		} else {
			read.csv(input$sp_file$datapath)
		}
	})

	sp_fit <- eventReactive(input$sp_fit, {
		dat <- sp_data()
		init <- c(input$sp_kd, input$sp_ka, input$sp_L0, input$sp_D0)
		res <- tryCatch(
			optim(par = init, fn = sp_error,
				  data = dat, DO_sat = input$sp_DO_sat,
				  method = "Nelder-Mead",
				  control = list(maxit = 2000),
				  hessian = TRUE),
			error = function(e) NULL
		)
		req(res)
		# Variance estimate from Hessian:
		# SSR-based: vcov ~= 2 * sigma^2 * solve(H), sigma^2 = SSR_min / (n - p)
		n_obs <- nrow(dat) * 2  # BOD + DO contribute
		p     <- length(res$par)
		sigma2 <- res$value / max(n_obs - p, 1)
		vcov <- tryCatch(2 * sigma2 * solve(res$hessian), error = function(e) NULL)
		ses  <- if (!is.null(vcov)) suppressWarnings(sqrt(diag(vcov))) else rep(NA, p)
		list(par = res$par, ses = ses, conv = res$convergence)
	})

	output$spRawPlot <- renderPlot({
		dat <- sp_data()
		ggplot(dat) +
			geom_point(aes(x = time, y = BOD), color = "red",       size = 2) +
			geom_point(aes(x = time, y = DO),  color = "steelblue", size = 2) +
			geom_hline(yintercept = input$sp_DO_sat, linetype = "dashed", color = "grey60") +
			labs(title = "Imported data (no model)",
				 x = "Time (days)", y = "Concentration (mg/L)") +
			theme_minimal(base_size = 13)
	})

	output$spFitText <- renderPrint({
		fit <- sp_fit()
		par <- fit$par
		ses <- fit$ses
		ci_kd <- c(par[1] - 1.96 * ses[1], par[1] + 1.96 * ses[1])
		ci_ka <- c(par[2] - 1.96 * ses[2], par[2] + 1.96 * ses[2])
		cat("Convergence code:", fit$conv, "(0 = success)\n\n")
		cat(sprintf("k_d = %.4f   95%% CI: [%.4f, %.4f]\n", par[1], ci_kd[1], ci_kd[2]))
		cat(sprintf("k_a = %.4f   95%% CI: [%.4f, %.4f]\n", par[2], ci_ka[1], ci_ka[2]))
		cat(sprintf("L0  = %.4f\n", par[3]))
		cat(sprintf("D0  = %.4f\n", par[4]))
	})

	output$spFitPlot <- renderPlot({
		fit <- sp_fit()
		dat <- sp_data()
		par <- fit$par
		fitted <- sp_simulate(L0 = par[3], D0 = par[4],
							  k_d = par[1], k_a = par[2],
							  times = seq(min(dat$time), max(dat$time), length.out = 200))
		fitted$DO <- input$sp_DO_sat - fitted$D
		ggplot() +
			geom_point(data = dat, aes(x = time, y = BOD), color = "red",       size = 2, alpha = 0.6) +
			geom_point(data = dat, aes(x = time, y = DO),  color = "steelblue", size = 2, alpha = 0.6) +
			geom_line (data = fitted, aes(x = time, y = L),  color = "red",       linewidth = 1) +
			geom_line (data = fitted, aes(x = time, y = DO), color = "steelblue", linewidth = 1) +
			labs(title = "Optimized Streeter-Phelps fit",
				 subtitle = "Points: data; Lines: fitted model",
				 x = "Time (days)", y = "Concentration (mg/L)") +
			theme_minimal(base_size = 13)
	})

	# ---------- Mass Transport 1D -------------------------------------------
	mt_sim <- eventReactive(input$mt_run, {
		dx <- 100
		L  <- 10000
		num_boxes <- L / dx + 1
		x_locations <- seq(0, L, by = dx)
		state <- numeric(num_boxes)
		state[max(1, min(num_boxes, input$mt_spike_pos))] <- input$mt_spike_C
		params_list <- list(u = input$mt_u, D = input$mt_D, dx = dx,
							C_upstream = input$mt_C_up, C_downstream = input$mt_C_down,
							bc_down = input$mt_bc)
		times <- seq(0, input$mt_tmax, by = 1)
		out <- ode(y = state, times = times, func = transport_model, parms = params_list)
		list(out = as.data.frame(out), x = x_locations)
	}, ignoreNULL = FALSE)

	output$mtHeatPlot <- renderPlot({
		sim <- mt_sim()
		df  <- sim$out
		x   <- sim$x
		mat <- as.matrix(df[, -1])
		rownames(mat) <- round(df$time, 2)
		colnames(mat) <- round(x, 0)
		mlt <- melt(mat)
		names(mlt) <- c("time", "distance", "Concentration")
		mlt$time     <- as.numeric(as.character(mlt$time))
		mlt$distance <- as.numeric(as.character(mlt$distance))
		ggplot(mlt, aes(x = distance, y = time, fill = Concentration)) +
			geom_tile() +
			scale_fill_viridis_c() +
			labs(title = "1D Mass Transport: Space-Time Diagram",
				 x = "Distance downstream (m)", y = "Time (days)") +
			theme_minimal(base_size = 13)
	})
}

shinyApp(ui, server, options = list(launch.browser = TRUE))
