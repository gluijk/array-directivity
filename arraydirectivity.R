# Calculating the Array factor (AF) and Directivity of linear Broadside arrays
# www.overfitting.net
# https://www.overfitting.net/2025/11/simulando-arrays-con-r.html


library(Cairo)

# Compute and plot the radiation pattern of a linear array
# of N identical radiating elements equally spaced
#
# This version ensures a true BROADSIDE definition (maximum perpendicular to the array axis)
#
# Conventions used:
# - Elements lie along the y-axis, centred at the origin.
# - 'theta' is the polar angle measured from broadside (= direction perpendicular to the array)
#   so theta = 0 points to the broadside direction. Positive theta rotates toward +y
# - Array Factor (AF) uses the phase term exp(i*k*x*cos(theta)) so that AF is maximum at theta=0
#
# Important features:
#  - AF is NOT normalized (you get absolute AF magnitude and AF dB = 20*log10(|AF|))
#  - Power directivity is computed as D(theta) = |AF(theta)|^2 / mean_theta(|AF|^2) and
#    returned both linear and in dB: 10*log10(D)
#  - The polar plot is drawn so that theta=0 points right (broadside) and increases clockwise
#    in degrees
#
# Args:
#   N: number of elements (integer)
#   d: element spacing in meters (numeric)
#   freq: frequency in Hz (numeric)
#   c: speed of sound in m/s (default 343)
#   theta.res: angular resolution in degrees (default 0.5)
#   maxdB: minimum dB shown in AF plot (for clipping visualization)
#
# Returns: invisible list with AF (complex), AF_mag, AF_power, directivity (linear & dB), angles, and x positions.


plot_array_pattern <- function(N = 2,
                               d = 0.30,
                               freq = 100,
                               c = 343,
                               theta.res = 0.1,
                               maxdB = -40,
                               name = "arrayradiationpattern.png") {
    
    if (N < 1) stop("N must be >= 1")
    if (d <= 0) stop("d must be > 0")
    if (freq <= 0) stop("freq must be > 0")
    
    # wavenumber
    k <- 2 * pi * freq / c
    
    # element positions along y-axis (centred at origin)
    idx <- seq(0, N - 1)
    y <- (idx - (N - 1)/2) * d
    elem_y <- y / (max(1, d*(N-1))) * 4  # plotting scaling of array
    
    # angles (theta=0 is broadside)
    theta_deg <- seq(-180, 180, by = theta.res)
    theta <- theta_deg * pi / 180
    
    # Array Factor (broadside maximum at theta=0)
    AF <- sapply(theta, function(t) sum(exp(1i * k * y * sin(t))))
    
    # magnitude and power
    AF_mag <- Mod(AF)
    AF_power <- AF_mag^2
    
    # AF magnitude in dB
    # (by adding eps we ensure the argument of log10() is never zero avoiding -Inf or errors)
    eps <- .Machine$double.eps

    # Directivity (linear)
    mean_power <- mean(AF_power)
    if (mean_power <= 0) mean_power <- eps
    directivity <- AF_power / mean_power
    
    # Polar coordinates for plotting
    r_af <- AF_mag
    x_af_plot <- r_af * cos(theta)   # horizontal = broadside
    y_af_plot <- r_af * sin(theta)   # vertical = array
    
    r_dir <- directivity           # for polar plot of linear directivity
    x_dir_plot <- r_dir * cos(theta)
    y_dir_plot <- r_dir * sin(theta)
    
    CairoPNG(name, width = 1280, height = 1080)
    
        # Plotting layout
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par))
        par(mfrow = c(2,2), mar = c(4,3,4,3), oma = c(3,2,3,2))
    
        # 1) AF magnitude (linear)
        plot(x_af_plot, y_af_plot, type='l', lwd=1, col='red', xlab='', ylab='',
             asp=1, main='AF magnitude (linear) = |AF|', cex.main=2, cex.axis=1.5,
             #xlim=c(-2,2), ylim=c(-2,2)
             )
    
        points(rep(0,length(y)), elem_y, pch=19, col='blue', cex=1.5)
        text(rep(0.05,length(y)), elem_y, labels=order(-seq_len(N)), pos=4, cex=1.5, col='blue')
        lines(c(0,0), c(elem_y[1],elem_y[N]), lwd=2, col='blue')
        text(min(x_af_plot), max(y_af_plot), labels=paste0("max |AF|=",round(max(AF_mag),1)),
             adj = c(0, 0), cex=2)
        abline(h=0, v=0)
        
        # 2) AF magnitude (dB)
        AF_mag_db <- 20 * log10(AF_mag + eps)
        plot(theta_deg, AF_mag_db, type='l', lwd=1, col='red',
             xlab='Angle (deg, 0 = broadside)', ylab='',
             main='AF magnitude (dB) = 20*log10(|AF|)',
             ylim=c(maxdB, max(AF_mag_db)+6), cex.main=2, cex.axis=1.5, cex.lab=1.2)
             #ylim=c(maxdB, 10), cex.main=2, cex.axis=1.5, cex.lab=1.2)
        text(min(theta_deg), max(AF_mag_db), labels=paste0("max |AF|(dB)=",round(max(AF_mag_db),1),"dB"),
             adj = c(0, 0), cex=2)
        abline(h=c(-12,-6,6,12), col='gray', lty='dotted')
        abline(h=0, col='gray')
        abline(v=0)
        
        # 3) Directivity (linear)
        plot(x_dir_plot, y_dir_plot, type='l', lwd=1, col='red', xlab='', ylab='', asp=1,
             main='Directivity (linear) = |AF|^2 / mean(|AF|^2)', cex.main=2, cex.axis=1.5,
             #xlim=c(-3.5,3.5), ylim=c(-2,2)
             )
        points(rep(0,length(y)), elem_y, pch=19, col='blue', cex=1.5)
        text(rep(0.05,length(y)), elem_y, labels=order(-seq_len(N)), pos=4, cex=1.5, col='blue')
        lines(c(0,0), c(elem_y[1],elem_y[N]), lwd=2, col='blue')
        text(min(x_dir_plot), max(y_dir_plot), labels=paste0("max D=",round(max(directivity),1)),
             adj = c(0, 0), cex=2)
        abline(h=0, v=0)
        
        # 4) Directivity (dB)
        directivity_db <- 10*log10(directivity + eps)
        plot(theta_deg, directivity_db, type='l', lwd=1, col='red',
             xlab='Angle (deg, 0 = broadside)', ylab='',
             main='Directivity (dB) = 10*log10(|AF|^2 / mean(|AF|^2))',
             ylim=c(maxdB, max(directivity_db)+6), cex.main=2, cex.axis=1.5, cex.lab=1.2)
        #ylim=c(maxdB, 10), cex.main=2, cex.axis=1.5, cex.lab=1.2)
        text(min(theta_deg), max(directivity_db), labels=paste0("max D(dB)=",round(max(directivity_db),1),"dB"),
             adj = c(0, 0), cex=2)
        abline(h=c(-12,-6,6,12), col='gray', lty='dotted')
        abline(h=0, col='gray')
        abline(v=0)
        
        # General title
        title(paste0('BROADSIDE vertical array: N=', N, ', d=', d, 'm, f=', freq, 'Hz'),
              outer = TRUE, cex.main = 2.5)
    
    dev.off()
        
    invisible(list(
        N=N, d=d, freq=freq, c=c, k=k, y=y,
        theta_deg=theta_deg, theta=theta, AF=AF, AF_mag=AF_mag,
        AF_mag_db=AF_mag_db, AF_power=AF_power,
        directivity=directivity, directivity_db=directivity_db
    ))
}


# Example 1: church loudspeakers
res <- plot_array_pattern(N = 5, d = 0.2, freq = 500)


# Example 2: d'Appolito configuration
for (f in seq(100, 3000, by=10)) {
    name=paste0("dappolito_", ifelse(f<1000,"0",""), f, "Hz.png")
    res <- plot_array_pattern(N=2, d=0.25, freq=f, name=name)
}

# First null frequency, single lobe array
f=343/(2*0.25)  # freq = v / (N * d)
res <- plot_array_pattern(N = 2, d = 0.25, freq = f)
