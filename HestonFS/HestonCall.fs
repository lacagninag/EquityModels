namespace HestonFS

open Fairmat.Math;
open Heston
open HestonEstimator
open DVPLI
open DVPLDOM
open System


module HestonContext = 

    type HestonInstance = {
        K : float
        Kappa : float
        Theta : float
        Sigma : float
        Rho : float
        S0 : float
        T : float
        V0 : float
        R : float
        Q : float
    }


module HestonCall =

    
    type HestonCall() = 

        static let phi (u: Complex) (x: HestonContext.HestonInstance) =
        
            let I = Complex.I
            let ss = x.Sigma * x.Sigma
            let tmp1 = I * x.Rho * x.Sigma * u
            let d = Complex.Sqrt((tmp1 - x.Kappa) ** 2.0 + ss * (I * u + u * u))
            let tmp2 = x.Kappa - tmp1
            let par = tmp2 - d
            let g = par / (tmp2 + d)

            let edT = Complex.Exp(-d * x.T)
            let numArg = Complex(1.0, 0.0) - g * edT
            let A = (x.Theta * x.Kappa) * (par * x.T - 2.0 * Complex.Log(numArg / (Complex(1.0, 0.0) - g))) / ss
            let B = x.V0 * (par * (Complex(1.0, 0.0) - edT) / numArg) / ss

            let log_s0 = Complex(Math.Log(x.S0), 0.0)
            let log_val = Math.Log(x.S0) +  x.R * x.T

            Complex.Exp(I * u * log_val + A + B)



        static let integrandFunc1 (u:double) (x: HestonContext.HestonInstance) = 
        
            let i = Complex.I
            let iu = i * u
            let A = Complex.Exp(-iu * Math.Log(x.K))

            (A * (phi (u - i) x) / iu ).Re

        static let integrandFunc2 (u:double) (x: HestonContext.HestonInstance) = 
        
            let i = Complex.I
            let iu = i * u
            let A = Complex.Exp(-iu * Math.Log(x.K))
            let uC = Complex(u, 0.0)

            (A * (phi uC x) / iu ).Re

        static let integrandFunc (u:double) (x: HestonContext.HestonInstance) = 
            integrandFunc1 u x - x.K * (integrandFunc2 u x)


        static let performIntegral (a: double) (b: double) (functionToIntegrate: TAEDelegateFunction1D) : double =
            
            let rec integrate x0 sum0 dt0 =
                let x = x0 + dt0
                let y0 = functionToIntegrate.Invoke x0
                let y1 = functionToIntegrate.Invoke x
                let sum = sum0 + 0.5 * (y0 + y1) * dt0
                let dt =
                    if x >= b then
                        b - x0
                    else
                        dt0 * 1.05
                if x >= b then
                    sum
                else
                    integrate x sum dt

            integrate a 0.0 (a / 10.0)

        

        static member HestonCallPrice (x:HestonContext.HestonInstance)=

            if Engine.Verbose > 0 then
                printfn "Pricing a call option with Heston model"
                printfn "Heston Parameters"
                printfn "kappa\ttheta\tsigma\trho\tv0"
                printfn "%f\t%f\t%f\t%f\t%f" x.Kappa x.Theta x.Sigma x.Rho x.V0
                printfn "Call Option Information"
                printfn "s0\tK\tT\tr\tq"
                printfn "%f\t%f\t%f\t%f\t%f" x.S0 x.K x.T x.R x.Q

            let F = x.S0 * Math.Exp((x.R - x.Q) * x.T)
            let firstTerm = 0.5 * (F - x.K)

            let a = 1E-8
            let b = 1000.0
            let functionToIntegrate = new TAEDelegateFunction1D ( fun u -> integrandFunc u x)
            
            let part1 = performIntegral a b functionToIntegrate
            let integral = part1 + a * functionToIntegrate.Invoke (a / 2.0)
            let unDiscountedCall = firstTerm + integral / Math.PI
            let call = Math.Exp(-x.R * x.T) * unDiscountedCall

            if Engine.Verbose > 0 then
                printfn "Call Price: %f" call

            call


