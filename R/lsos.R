lsos <-
function (..., n = 10) 
{
    .ls.objects(..., order.by = "Size", decreasing = TRUE, head = TRUE, 
        n = n)
}
