namespace HestonFS


open Fairmat.Math;
open Heston
open HestonEstimator
open DVPLI
open DVPLDOM
open System
open HestonFS.HestonContext




module HestonGreekF =

    type HestonNumericalGreeks() =

        static let GreeksBumper (bumpPercentage: float) (init: float) (f: Func<float, float>) (unbumpedPrice: float option) =
            let deltaX = bumpPercentage * init

            match unbumpedPrice with
            | Some price ->
                let bumpedPrice = f.Invoke(init + deltaX)
                (bumpedPrice - price) / deltaX
            | None ->
                let bumpedPrice = f.Invoke(init + deltaX)
                let unbumpedPrice = f.Invoke(init - deltaX)
                (bumpedPrice - unbumpedPrice) / (2.0 * deltaX) 


        static member VegaCall (x:HestonInstance) (bumpPercentage : double)= 

            let verbosity = Engine.Verbose

            if verbosity > 0 then
                printfn "Calculating vega with Heston model"
                printfn "Heston Parameters"
                printfn "kappa\ttheta\tsigma\trho\tv0"
                printfn "%f\t%f\t%f\t%f\t%f" x.Kappa x.Theta x.Sigma x.Rho x.V0
                printfn "Call Option Information"
                printfn "s0\tK\tT\tr\tq"
                printfn "%f\t%f\t%f\t%f\t%f" x.S0 x.K x.T x.R x.Q

            let callPrice (initialVariance: double) = 
                HestonCall.HestonCallPrice(kappa = x.Kappa, theta = x.Theta, sigma = x.Sigma, rho = x.Rho, 
                                            v0 = initialVariance, s0 = x.S0, T = x.T, K = x.K, r = x.R, q = x.Q)

                                
            Engine.Verbose <- 0
            let vega = 
                GreeksBumper 
                    bumpPercentage 
                    x.V0 
                    callPrice
                    None
            Engine.Verbose <- verbosity

            if verbosity > 0 then
                printfn "Vega"
                printfn "%f" vega

            vega


        static member DeltaGammaCall  (x: HestonInstance) (bumpPercentage: double option) (unBumpedPrice : double option) : (double*double) = 
            
            let price = 
                match unBumpedPrice with
                | Some p -> p
                | None -> HestonCall.HestonCallPrice(kappa = x.Kappa, theta = x.Theta, sigma = x.Sigma, rho = x.Rho, v0 = x.V0, s0 = x.S0, T = x.T, K = x.K, r=x.R, q=x.Q)

            let bump =
                match bumpPercentage with
                | Some b -> b
                | None -> 0.01

            let deltaS = bump * x.S0

            let callPriceBumpUp = HestonCall.HestonCallPrice(kappa = x.Kappa, theta = x.Theta, sigma = x.Sigma, rho = x.Rho, v0 = x.V0, s0 = x.S0 + deltaS, T = x.T, K = x.K, r=x.R, q=x.Q)
            let callPriceBumpDown = HestonCall.HestonCallPrice(kappa = x.Kappa, theta = x.Theta, sigma = x.Sigma, rho = x.Rho, v0 = x.V0, s0 = x.S0 - deltaS, T = x.T, K = x.K, r=x.R, q=x.Q)



            (callPriceBumpUp - callPriceBumpDown) / (2.0 * deltaS), (callPriceBumpUp - 2.0 * price + callPriceBumpDown) / (deltaS * deltaS)

        



        

        

